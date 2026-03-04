from pipeline import PipelineTemplate
from typing import Callable
from molscrub import Scrub
from tqdm import tqdm
import pandas as pd
from functools import partial
import multiprocessing as mp
from rdkit import Chem
from meeko import PDBQTReceptor
from meeko import PDBQTMolecule
from meeko import ResidueChemTemplates
from meeko import Polymer
from meeko import MoleculePreparation
from pathlib import Path
from os import path,makedirs
from prody import parsePDB
from meeko import PDBQTWriterLegacy
from meeko.analysis import FingerprintInteractions
from lib import *
import prody
from random import randint
from collections import defaultdict
import json
from time import time
import os



prody.confProDy(verbosity="none")

def calc_interactions(ligand_id:str,ligand_file:str,target_id:str,target_file:str):
    ligand = PDBQTMolecule.from_file(ligand_file)
    target = PDBQTReceptor.from_pdbqt_filename(target_file)

    calculator = FingerprintInteractions(target)

    calculator.run(ligand)

    interactions = calculator.to_dataframe()
    
    interaction_types = [i for i, _ in interactions.columns]
    residues = [r for _, r in interactions.columns]

   
    df_t = pd.DataFrame({
        'InteractionType': interaction_types,
        'Residue': residues,
        'Value': interactions.iloc[0].values 
    })

    interactions_filtered = df_t[df_t['Value'] == 1]
    
    agg = interactions_filtered.groupby('InteractionType')['Residue'].agg(list)

    
    result_df = pd.DataFrame([agg])


    result_df["ligand_id"] = ligand_id
    result_df["target_id"] = target_id
    result_df["file"] = ligand_file.replace("pdbqt","sdf")

    result_dict = result_df.to_dict(orient='records')[0]
    
    return result_dict
    
    


def prepare_ligand(mol,ligand_id:str,target_id:str,n_conf:int,out_dir:str,align_pivot:Iterable,hydrate:bool=False,save_sdf:bool=True,preparation_config:dict=None,*args,**kwargs):
    mol = Chem.Mol(mol)

    ligand_dir=path.join(out_dir,ligand_id)
    
    # Molecule preparation from Meeko
    if preparation_config is not None:
        preparator = MoleculePreparation.from_config(preparation_config) #load config
    else:   
        preparator = MoleculePreparation(hydrate=hydrate) #hydrate protocol

    makedirs(ligand_dir,exist_ok=True)

    Chem.rdMolAlign.AlignMolConformers(mol)
    
    

    sizes=[]
    confs=[]
    sdf=[]
    for conf_id in range(mol.GetNumConformers()):
        
        center,conf_size=get_box(mol=mol,conf_id=conf_id)
        aligned=coords_align(mol,align_pivot,center,conf_id)
        
        sizes.append(conf_size)
    
        
         
        conf=preparator.prepare(aligned,conformer_id=conf_id)
        
        conf_pdbqt = PDBQTWriterLegacy.write_string(conf[0])
        write_pdbqt(dir=ligand_dir,filename=f"{ligand_id}_{conf_id}",pdbqt=conf_pdbqt) #save conf pdbqt
        
        
        
        
        confs.append(f"{ligand_dir}/{ligand_id}_{conf_id}.pdbqt")
        if save_sdf:
            sdf_file = path.join(ligand_dir,f"{ligand_id}_{conf_id}.sdf")
            sdf.append(sdf_file)
            writer = Chem.SDWriter(sdf_file)
            writer.write(aligned, confId=conf_id)
            

    writer.close()
    sizes=np.array(sizes)
    opt_size=sizes.max(axis=0)

    if save_sdf:
        return {"ligand_id":f"{ligand_id}","target_id":target_id,"conf_pdbqt":tuple(confs),"conf_sdf":tuple(sdf),"box_center":align_pivot,"box_size":opt_size}
    else:
        return {"ligand_id":f"{ligand_id}","target_id":target_id,"conf_pdbqt":tuple(confs),"box_center":align_pivot,"box_size":opt_size}

    


def prepare_target(file:str,
    target_id:str,
    altloc:str,
    box:dict,
    hydrate:bool=False,
    charge_type="gesteiger",
    config:dict=None,
    backend:str="prody",
    *args,**kwargs):

    
    outdir = str(Path(file).parent)
    outpath=path.join(outdir,f"{target_id}.pdbqt")
    box={"box_center":box[0],"box_size":box[1]}
    templates = ResidueChemTemplates.create_from_defaults() #create from defaults for now
    if config is not None:
        mk_prep = MoleculePreparation.from_config(config)
    else:
        mk_prep = MoleculePreparation(hydrate=hydrate,charge_model=charge_type)

    #load pdb with H
    if backend.lower() == "prody":
        target=parsePDB(file)

        polymer=Polymer.from_prody(
        target,
        templates,  # residue_templates, padders, ambiguous,
        mk_prep,
        #set_template,
        #delete_residues,
        allow_bad_res=True,
        #blunt_ends=True,
        #wanted_altloc=wanted_altloc,
        default_altloc=altloc
    )

    elif backend.lower() == "file":
        with open(file,"r") as pdb_file:
            pdb_string=pdb_file.read()
        
        polymer=Polymer.from_pdb_string(
        pdb_string,
        templates,  # residue_templates, padders, ambiguous,
        mk_prep,
        #set_template,
        #delete_residues,
        allow_bad_res=True,
        #blunt_ends=blunt_ends,
        #wanted_altloc=wanted_altloc,
        default_altloc=altloc,
    )
    else:
        raise ValueError(f"Invalid backend:{backend}")

    pdbqt_tuple = PDBQTWriterLegacy.write_from_polymer(polymer)
    
    pdbqt=write_pdbqt(dir=outdir,filename=f"{target_id}",pdbqt=pdbqt_tuple)

    with open(path.join(outdir,f"{target_id}.json"),"w") as json:
        json.write(polymer.to_json())

    
    return {"target_id":target_id,"path":outpath,"pdb_file":file,"json":path.join(outdir,f"{target_id}.json")} | box #return dict with target and box specs
    


class MeekoPipeline(PipelineTemplate):
   
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
                        "hydrate": hydrate,
                        "output": output,
                        "save_sdf": save_sdf,
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
        hydrate:bool=False,
        charge_type="gasteiger",
        config:dict=None,
        backend:str="prody",
        box:dict=None,
        chunk_size:int=None,
        *args,
        **kwargs
        ):

        targets = self.__getattr__("target_data")   
        #targets.set_index("target_id", inplace=True)

        
        tasks = (
            (row.file, row.target_id, row.altloc, box[row.target_id])
            for row in targets.itertuples(index=False)
        )

        prep_f = partial(prepare_target, hydrate=hydrate, charge_type=charge_type, config=config, backend=backend)

        prepared_sets = []
        n_workers = min(self.n_proc, targets.shape[0]) 

        with mp.Pool(n_workers) as pool:
            for r in tqdm(pool.imap_unordered(partial(map_wrapper, func=prep_f), tasks),
                        total=targets.shape[0],
                        desc="Target preparation",
                        dynamic_ncols=True):
                
                prepared_sets.append(r)

        prepared = pd.merge(pd.DataFrame(prepared_sets),targets, on="target_id", how="left")
        self.__setattr__("target_data", prepared.set_index("target_id"))
        print(f"\nSuccessfully prepared {prepared.shape[0]} target(s)\n")
        


    def init_docking(self, safety_factor: float = 0.9, **docking_kwargs):

        usable_cores = max(1, int(self.n_proc * safety_factor))


        self.__setattr__("docking_n_work",int(usable_cores/2))
        docking_kwargs["threads"] = 2
        """
        if self.ligands.shape[0] > usable_cores:
            self.__setattr__("docking_n_work",int(usable_cores/2))
            docking_kwargs["cpu"] = 2

        else:
            self.__setattr__("docking_n_work",self.ligands.shape[0])
            docking_kwargs["cpu"] = max(1,usable_cores//self.ligands.shape[0])

        """
        self.__setattr__("docking_config",dict(docking_kwargs))
        self.__setattr__("docking_init",True)


 
        
    
    def calculate_maps(self,map_f:Callable,map_mode:str,target_id:str,fit_ligands:bool=False,target_box:bool=False,**map_kwargs):
        
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
        
            
            print(f"Calculating maps for {target_id}")
            coupled=map_f(target=target_path,
            target_id=target_id,
            box_center=box_center,
            box_size = box_size,
            save_prefix=save_prefix,
            **self.docking_config)

            map_sets=[]
            for _id in tqdm(ligands["ligand_id"].tolist(),desc=f"Generating tasks"):
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
        pass
        

    def screen_single_target(self,dock_f:Callable,target_id:str,report_fn:str=None,**dock_kwargs):

        
        screen_path = f"{self.result_workpath}/{target_id}"
        makedirs(screen_path,exist_ok=True)

        if not report_fn:
            report_fn = f"{screen_path}/{target_id}_{self.run_name}"
        else:
            report_fn = f"{screen_path}/{report_fn}"

        self.__setattr__("docking_report",report_fn)

        n_ligands = self.ligand_data.shape[0]
        
        
        
        total = pd.merge(
            self.ligand_data.drop(columns=["box_center", "box_size"]),
            self.map_data,
            on=["ligand_id", "target_id"],
            how="left"
        )

        
        with tqdm(total=self.replicates,desc=f"Docking {target_id}, replicates {self.replicates}",position=0) as replicate_bar:

            

            dock_sets = []
   
            for replicate in range(1, self.replicates + 1):
               

                replicate_path = f"{screen_path}/replicate {replicate}"
                makedirs(replicate_path, exist_ok=True)

                if not dock_kwargs.get("seed"):
                    dock_kwargs["seed"] = randint(100000, 999999)

                self.__setattr__("docking_config", dock_kwargs)

                dock_func = partial(dock_f, **self.docking_config)

                
                    
                tasks = (
                    (
                        row.ligand_id,
                        conf_path,
                        row.target_id,
                        row.coupled_fn,
                        replicate_path,
                        replicate,
                        conf_id
                    )
                    for row in total.itertuples(index=False)
                    for conf_id, conf_path in enumerate(row.conf_pdbqt)
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

    def analyze_interactions(self,target_id:str):
        target_path = str(self.targets.loc[target_id,"path"])
        interaction_calculator = partial(calc_interactions,target_id=target_id,target_file=target_path)

        interaction_sets=[]
        

        tasks = ((row.ligand_id, row.file.replace("sdf", "pdbqt")) 
                for row in self.report.itertuples(index=False))

        n_workers = min(self.n_proc, self.report.shape[0]) if tasks else 1

        with mp.Pool(n_workers) as pool:
            for result in tqdm(pool.imap_unordered(partial(map_wrapper,func=interaction_calculator), tasks),total=self.report.shape[0],desc=f"Calculating interactions for {target_id}"):
                interaction_sets.append(result)

        analysis = pd.DataFrame(interaction_sets)
        analysis = analysis.drop_duplicates(
        subset=["ligand_id", "target_id", "file"]
        )
        analysis_report = pd.merge(self.report,analysis,on=["ligand_id","target_id","file"],how="left")
        self.__setattr__("report",analysis_report)
            
        analysis_report.to_csv(f"{self.result_workpath}/{target_id}/{target_id}_{self.run_name}_interactions.csv",encoding="utf-8",sep=",",index=False)
                
    @property
    def create_report(self):
        print("Creating report")
        interaction_cols = ["HBAcceptor", "Hydrophobic", "HBDonor"]

        output_file = path.join(self.result_workpath,f"{self.exp_id}_{self.run_name}.xlsx")
        interaction_cols = ["HBAcceptor", "Hydrophobic", "HBDonor"]

        report = self.__getattr__("full_report")
        report.drop(columns=['file'], inplace=True)
        main_cols = [
            "ligand_id",
            "target_id",
            "energy (kcal/mol)",
            "temperature (Celsius)",
            "Ki (microM)",
            "Ligand Eff",
            "reproduced",
            "max RMSD",
            "redocked",
        ]

        with pd.ExcelWriter(output_file,engine="xlsxwriter") as writer:

            for target_id, target_df in report.groupby("target_id"):

               
                target_df = target_df.copy()

                
                for col in interaction_cols:
                    if col in target_df.columns:
                        target_df[col] = target_df[col].apply(
                            lambda val: (
                                ", ".join(str(x) for x in val)
                                if isinstance(val, list)
                                else val.strip("[]").replace("'", "")
                                if isinstance(val, str)
                                else val
                            )
                        )

                
                main_df = target_df[main_cols].copy()

               
                interaction_df = target_df[
                    ["ligand_id", "target_id", "energy (kcal/mol)"] + interaction_cols
                ].copy()

                main_sheet = f"{target_id}_summary"[:31]
                inter_sheet = f"{target_id}_interactions"[:31]

                main_df.to_excel(writer, sheet_name=main_sheet, index=False)
                interaction_df.to_excel(writer, sheet_name=inter_sheet, index=False)



        
            





    

