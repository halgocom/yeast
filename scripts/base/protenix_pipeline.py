from pipeline import PipelineTemplate
from typing import Callable
from molscrub import Scrub
from tqdm import tqdm
import pandas as pd
from functools import partial
import multiprocessing as mp
from rdkit import Chem
from pathlib import Path
from os import path,makedirs
from prody import parsePDB
from lib import *
import prody
from random import randint
from collections import defaultdict
import pdb2pqr




prody.confProDy(verbosity="none")
def prepare_target(file:str,
    target_id:str,
    box:dict,
    *args,**kwargs):

    box={"box_center":box[0],"box_size":box[1]}


    
    return {"target_id":target_id,"path":file} | box #return dict with target and box specs



def prepare_ligand(mol,ligand_id:str,target_id:str,n_conf:int,out_dir:str,align_pivot:Iterable,hydrate:bool=False,save_sdf:bool=True,preparation_config:dict=None,*args,**kwargs):
    mol = Chem.Mol(mol)

    ligand_dir=path.join(out_dir,ligand_id)


    makedirs(ligand_dir,exist_ok=True)

    Chem.rdMolAlign.AlignMolConformers(mol)
    
    

    sizes=[]
    confs=[]
    sdf=[]
    for conf_id in range(mol.GetNumConformers()):
        
        center,conf_size=get_box(mol=mol,conf_id=conf_id)
        aligned=coords_align(mol,align_pivot,center,conf_id)
        
        sizes.append(conf_size)
    
        
         

        sdf_file = path.join(ligand_dir,f"{ligand_id}_{conf_id}.sdf")
        sdf.append(sdf_file)
        writer = Chem.SDWriter(sdf_file)
        writer.write(aligned, confId=conf_id)
            

    writer.close()
    sizes=np.array(sizes)
    opt_size=sizes.max(axis=0)

    
    return {"ligand_id":f"{ligand_id}","target_id":target_id,"conf_sdf":tuple(sdf),"box_center":align_pivot,"box_size":opt_size}


class ProtenixPipeline(PipelineTemplate):
    
    def optimize_ligand(self,opt_f:Callable=None,n_conf:int=1,protonate:bool=True,pH_min:float=None,pH_max:float=None,*opt_args,**opt_kwargs):
        ligands = self.__getattr__("ligand_data")
        
        if n_conf < 1 or (not isinstance(n_conf,int)):
            raise ValueError(f"Unknown number of configurations to generate:{n_conf}")
        self.__setattr__("N_CONF_LIGAND",n_conf)


        opt_kwargs["n_conf"] = n_conf
        opt_f = partial(opt_f, *opt_args, **opt_kwargs)
        func_name = opt_f.func.__name__

        
        tasks = ((row.mol, row.id) for row in ligands.itertuples(index=False))

        optimized_sets = []
        n_workers = min(self.n_proc, ligands.shape[0])

        with mp.Pool(n_workers) as pool:
            for r in tqdm(pool.imap_unordered(partial(map_wrapper, func=opt_f), tasks),
                        total=ligands.shape[0],
                        desc=f"Optimizing geometry with {func_name}"):
                
                optimized_sets.append(r)


        
        if protonate:
            if (pH_min) and (not pH_max):
                pH_max = pH_min
            
            if (pH_max) and (not pH_min):
                pH_min = pH_max
            
            if (not pH_max) and (not pH_min):
                pH_max = 7.4
                pH_min = 7.4

            
            scrub=Scrub(ph_low=pH_min,ph_high=pH_max,skip_gen3d=True,debug=True)
                
            try:  
                if n_conf ==1:
                    scrubbed_mols=[{"id":ligand["id"],"mol":scrub(ligand["mol"])[0]} for ligand in tqdm(optimized_sets,desc=f"Protonation pH range:{pH_min}-{pH_max}",position=0)]
                else:
                    scrubbed_mols=[]
                    
                    for ligand in tqdm(optimized_sets,desc=f"Protonation pH range:{pH_min}-{pH_max}",position=0):
                        new_mol =Chem.Mol(ligand["mol"])
                        new_mol.RemoveAllConformers()
                        new_mol.SetProp("id",ligand["id"])
                        for conf_id in range(n_conf):
                            conf = Chem.Conformer(ligand["mol"].GetConformer(id=conf_id))
                            conf_mol = Chem.Mol(ligand["mol"])
                            conf_mol.RemoveAllConformers()
                            conf_mol.AddConformer(conf, assignId=True)
                            scrubbed_conf = scrub(conf_mol)[0]
                            new_mol.AddConformer(scrubbed_conf.GetConformer(),assignId=True)

                        #Chem.rdMolAlign.AlignMolConformers(new_mol)
                        scrubbed_mols.append({"id":ligand["id"],"smiles":Chem.MolToSmiles(new_mol),"mol":new_mol})
                
                ligands = pd.DataFrame(scrubbed_mols)
            except Exception:
                pass
        else:
            print("\nSkipping protonation\n")
            ligands = pd.DataFrame(optimized_sets)


        self.__setattr__("ligand_data",ligands) #join all results together
        print(f"\nSuccessfully optimized {ligands.shape[0]} ligands(s)\n")



    def process_ligands(self,hydrate:bool=False,output:str="file",save_sdf:bool=True,preparation_config:dict=None,*args,**kwargs):
        ligands = self.__getattr__("ligand_data")
        n_target = self.__getattr__("target_data").shape[0]
        n_ligands = ligands.shape[0]



    
        tasks = ((row.mol.ToBinary(), row.id) for row in ligands.itertuples(index=False))

        prepared_sets = []

        with tqdm(total=n_target, desc=f"Preparing ligands for {n_target} target(s)", position=0) as target_bar:
            for target_id in self.target_data.index:
                with tqdm(total=n_ligands, desc=f"Target {target_id}", position=1, leave=False, dynamic_ncols=True) as ligand_bar:
                    
                    target_box_center, target_box_size = self.target_data.loc[target_id, ["box_center", "box_size"]].to_numpy()

                    prep_kwargs = {
                        "n_conf": self.N_CONF_LIGAND,
                        "out_dir": self.ligands_workpath,
                        "output": output,
                        "preparation_config": preparation_config,
                        "align_pivot": target_box_center,
                        "target_id": target_id
                    }

                    prep_f = partial(prepare_ligand, **prep_kwargs)

                    n_workers = min(self.n_proc, ligands.shape[0]) 
                    with mp.Pool(n_workers) as pool:
                        for r in pool.imap_unordered(partial(map_wrapper, func=prep_f), tasks):
                            
                            prepared_sets.append(r)
                            ligand_bar.update(1)

                target_bar.update(1)

        prepared = pd.DataFrame(prepared_sets)
        print(f"Successfully prepared {prepared.shape[0]} ligand(s)")   
        self.__setattr__("ligand_data",prepared) #join all results together
        prepared.to_csv(f"{self.ligands_workpath}/{self.exp_id}_{self.run_name}",encoding="utf-8",sep=",",index=False)


    def process_targets(self,
        backend:str="prody",
        box:dict=None,
        chunk_size:int=None,
        *args,
        **kwargs
        ):

        targets = self.__getattr__("target_data")   
        #targets.set_index("target_id", inplace=True)

        
        tasks = (
            (row.file, row.target_id,box[row.target_id])
            for row in targets.itertuples(index=False)
        )

        

        prepared_sets = []
        n_workers = min(self.n_proc, targets.shape[0]) 

        with mp.Pool(n_workers) as pool:
            for r in tqdm(pool.imap_unordered(partial(map_wrapper, func=prepare_target), tasks),
                        total=targets.shape[0],
                        desc="Target preparation",
                        dynamic_ncols=True):
                
                prepared_sets.append(r)

        prepared = pd.merge(pd.DataFrame(prepared_sets),targets, on="target_id", how="left")
        self.__setattr__("target_data", prepared.set_index("target_id"))
        print(f"\nSuccessfully prepared {prepared.shape[0]} target(s)\n")


    def init_docking(self,safety_factor:float=0.9,**docking_kwargs):
        
        docking_kwargs.setdefault("threads", 1)
        usable_cores = max(1, int(self.n_proc * safety_factor))


        self.__setattr__("docking_n_work",int(usable_cores/2))
        docking_kwargs["threads"] = 2

        self.__setattr__("docking_config",docking_kwargs)
        self.__setattr__("docking_init",True)



    def set_maps(self,map_f:Callable,map_mode:str,target_id:str,fit_ligands:bool=False,target_box:bool=False,**map_kwargs):
        if not self.docking_init:
            raise RuntimeError("This error is a safeguard,use init_docking method to overule it")

        
        ligands= self.__getattr__("ligand_data")
        targets=self.__getattr__("target_data")
        
        n_ligands = ligands.shape[0]
        
        
        box={}
        save_prefix=f"{self.maps_workpath}/{target_id}"
        makedirs(save_prefix,exist_ok=True)
        target_box_center,target_box_size,target_path = targets.loc[target_id,["box_center","box_size","path"]].to_numpy()
        #try:
        if not map_mode:
            raise RuntimeError(f"Missing map mode for {target_id}")
        
        if map_mode == "o4a":
            if fit_ligands and target_box:
                raise RuntimeError("fit_ligands and target_box option can't both be used")
            
            if (not fit_ligands) and (not target_box):
                raise RuntimeError("At least one between fit_ligands and target_box must be choosen")
            
            
            if target_box:
                box_center=np.array(target_box_center)
                box_size=np.array(target_box_size)
            
            if fit_ligands:
                box_center=target_box_center
                box_dims=np.vstack(ligands.loc[ligands["target_id"].eq(target_id),"box_size"].to_numpy()) #get all box dims and convert complatiple array
                box_size=np.max(box_dims, axis=0)
        
            
            
            coupled=map_f(target=target_path,
            target_id=target_id,
            box_center=box_center,
            box_size = box_size,
            save_prefix=save_prefix,
            **self.docking_config)

            map_sets=[]
            for _id in tqdm(ligands["ligand_id"].tolist(),desc=f"Generating map tasks"):
                entry = coupled.copy()
                entry["ligand_id"] = _id
                map_sets.append(entry)
            
            
            
        if map_mode == "o4e":

           
            
            map_sets = []
            map_calculator = partial(map_f, **self.docking_config)

            
            tasks = (
                (
                    target_path,
                    target_id,
                    target_box_center,
                    row.box_size,
                    f"{save_prefix}/{row.ligand_id}",
                    row.ligand_id
                )
                for row in ligands.itertuples(index=False)
            )

            

            with mp.Pool(self.docking_n_work) as pool:
                for r in tqdm(pool.imap_unordered(partial(map_wrapper, func=map_calculator), tasks),total=ligands.shape[0],desc=f"Map calculation for {target_id}"):
                    map_sets.append(r)
    
            
    
        self.__setattr__("map_data",pd.DataFrame(map_sets))

    def screen_single_target(self,dock_f:Callable,target_id:str,report_fn:str=None,**dock_kwargs):

        
        screen_path = f"{self.result_workpath}/{target_id}"
        makedirs(screen_path,exist_ok=True)
        
        target_box_center,target_box_size,target_path = self.targets.loc[target_id,["box_center","box_size","path"]].to_numpy()
        
        if not report_fn:
            report_fn = f"{screen_path}/{target_id}_{self.run_name}"
        else:
            report_fn = f"{screen_path}/{report_fn}"

        self.__setattr__("docking_report",report_fn)

        n_ligands = self.ligand_data.shape[0]
        
        
        
        #total = pd.merge(
            #self.ligand_data.drop(columns=["box_center", "box_size"]),
            #self.map_data,
            #on=["ligand_id", "target_id"],
            #how="left"
        #)

        
        with tqdm(total=self.replicates,desc=f"Docking {target_id}, replicates {self.replicates}",position=0) as replicate_bar:

            

            dock_sets = []
            for replicate in range(1, self.replicates + 1):
               
                seed = randint(10000000,999999999999)
                replicate_path = f"{screen_path}/replicate {replicate}"
                makedirs(replicate_path, exist_ok=True)

                task_rep = replicate
                task_rep_path = replicate_path


                

                dock_func = partial(dock_f, **self.docking_config)

                #receptor_id:str,receptor:str,ligand:str,ligand_id:str,replicate_id:int,save_path:str,box_center:Iterable,box_size:Iterable,conf_id:int,threads:int)
                tasks = (
                    (
                        target_id,
                        target_path,
                        conf_path,
                        row.ligand_id,
                        task_rep,
                        task_rep_path,
                        target_box_center,
                        target_box_size,
                        conf_id,
                        seed
                    )
                    for row in self.ligands[:1].itertuples(index=False)
                    for conf_id, conf_path in enumerate(row.conf_sdf)
                )


                total_tasks = self.N_CONF_LIGAND*self.ligands.shape[0]

                with tqdm(total=total_tasks,position=1,leave=False,desc=f"Docking {target_id}, {self.ligands.shape[0]} ligands, replicate {replicate}") as ligand_bar:

                    with mp.Pool(self.docking_n_work) as pool:

                        for r in pool.imap_unordered(partial(map_wrapper, func=dock_func),tasks):

                            dock_sets.append(r)

                            ligand_bar.update(1)

                replicate_bar.update(1)

            grouped = defaultdict(list)


            for r in dock_sets:
                key = (r['ligand_id'], r['target_id'], r['replicate'])
                grouped[key].append(r)

            final_list = []
            for key, conf_list in grouped.items():
            # ordina solo se vuoi (qui lasciamo l'ordine originale)
                merged = {
                    'ligand_id': key[0],
                    'target_id': key[1],
                    'replicate': key[2],
                    'docked_pdbqt': [r['docked_pdbqt'] for r in conf_list],
                    'docked_sdf': [r['docked_sdf'] for r in conf_list],
                    'energies': [r['energies'] for r in conf_list]
                }
                final_list.append(merged)

                    
                    
            docked = pd.DataFrame(final_list)
                
            docked.to_csv(f"{report_fn}.csv", sep=",", index=False, encoding="utf-8")
            self.__setattr__("docking_data",docked)
    pass
        
    


   



    