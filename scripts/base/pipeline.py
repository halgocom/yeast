from abc import ABC,abstractmethod
from dataclasses import dataclass
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from prody import parsePDB,saveAtoms,writePDB,AtomGroup
from os import getcwd,listdir,path,mkdir,makedirs,getenv
from typing import Iterable
from moloader import *
import pandas as pd
from typing import Callable
import subprocess
import sys
import time
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
from lib import *
import ast
from dotenv import load_dotenv
from rdkit import RDLogger
from rdkit.Chem import MolFromSmiles,MolFromPDBFile,Mol,rdMolAlign,SDMolSupplier,SanitizeMol,AllChem
from collections import defaultdict
from functools import reduce
import operator
import sys
from posebusters import PoseBusters
from pathlib import Path
from math import exp
from pdb2pqr import main as pdb2pqr_main

 



def inibition_costant(bind_energy:float,temp:float,energy_unit:str,*args):
    gas_costants={"kcal":0.00198588,"kj":0.00831447}
    R = gas_costants.get(energy_unit)
    
    ki = exp((bind_energy/(R*temp)))
    return ki*1000000

def ligand_efficiency(bind_energy:float,ligand_file:str):
    mol = SDMolSupplier(ligand_file)[0]
    n_heavy = mol.GetNumHeavyAtoms()

    return (-1)*(bind_energy/n_heavy)


def calc_parameters(ligand_id:str,ligand_file:str,target_id:str,bind_energy:float,temp:float,energy_unit:str,round_digit:int,parameters:Iterable):

    temp_k = temp + 273.15

    par_log = {"ligand_id":ligand_id,"target_id":target_id,f"energy ({energy_unit}/mol)":bind_energy,"temperature (Celsius)":temp}
    #parameters_f = {"ki":inibition_costant,"le":ligand_efficiency}

    for par in parameters:

        if par == "ki":
            value = np.round(inibition_costant(bind_energy,temp_k,energy_unit),round_digit)
            par_log.update({"Ki (microM)":value})
        elif par == "le":
            value = np.round(ligand_efficiency(bind_energy,ligand_file),round_digit)
            par_log.update({"Ligand Eff":value})
        else:
            par_log.update({par:None})

    return par_log





def bust_poses(ligand_id:str,dock_results:list,dock_energies:list,replicate:int,target_id:str,target_file:str,strict_clash_check:bool=False):
    
    
    
    if strict_clash_check:
        check_columns=['mol_pred_loaded', 'mol_cond_loaded', 'sanitization',
        'inchi_convertible', 'all_atoms_connected', 'no_radicals',
        'bond_lengths', 'bond_angles', 'internal_steric_clash',
        'aromatic_ring_flatness', 'non-aromatic_ring_non-flatness',
        'double_bond_flatness', 'internal_energy',
        'protein-ligand_maximum_distance', 'minimum_distance_to_protein',
        'minimum_distance_to_organic_cofactors',
        'minimum_distance_to_inorganic_cofactors', 'minimum_distance_to_waters',
        'volume_overlap_with_protein', 'volume_overlap_with_organic_cofactors',
        'volume_overlap_with_inorganic_cofactors',
        'volume_overlap_with_waters']
    else:
        check_columns=[
        'mol_pred_loaded', 'mol_cond_loaded', 'sanitization',
        'inchi_convertible', 'all_atoms_connected', 'no_radicals',
        'bond_lengths', 'bond_angles', 'internal_steric_clash',
        'aromatic_ring_flatness', 'non-aromatic_ring_non-flatness',
        'double_bond_flatness', 'internal_energy',
        'protein-ligand_maximum_distance', 'minimum_distance_to_protein','volume_overlap_with_protein']

    clash_checker = PoseBusters(config="dock")
    
    target_mol = MolFromPDBFile(target_file,sanitize=False,removeHs=False)
    docked_clash = clash_checker.bust(mol_pred=dock_results,mol_cond=target_mol)
    docked_clash["energy"] = pd.Series(dock_energies, index=docked_clash.index)
    

    dock_mask = np.all(docked_clash[check_columns].to_numpy(dtype=bool), axis=1)
    
    best_confs = docked_clash.loc[dock_mask]#, ["file","energy"]]

    file = best_confs["energy"].idxmin()
    energy = best_confs.loc[file, "energy"]


    
    bust_log ={"replicate":replicate,
    "ligand_id":ligand_id,
    "target_id":target_id,
    "file":file[0],
    "energy":energy}

    return bust_log
    
def check_producibility(ligand_id:str,files:Iterable,energy:Iterable,target_id:str,threshold:float=2.0):
    status = "No"

    best_energy = min(energy)
    test_idx = energy.index(best_energy)
    
    best_file =files[test_idx]

    best_replicate_mol = SDMolSupplier(files[test_idx])[0]

    test_files = files[:test_idx] + files[test_idx+1:]

    test_mol = [SDMolSupplier(file)[0] for file in test_files]
    
    rmsd=[rdMolAlign.GetBestRMS(best_replicate_mol,mol) for mol in test_mol]

    if all(r <= threshold for r in rmsd):
        status = "Yes"

    return {"ligand_id":ligand_id,"file":best_file,"target_id":target_id,"energy":best_energy,"reproduced":status,"max RMSD":max(rmsd)}


    




def target_clean(filepath:str,target_id:str,altloc:str=None,chain:str=None,cryst_ligand:str=None,pH:float=7.4,forcefield:str=None,keepWater:bool=False,repair:bool=True,out_path:str=None):
    
    
    outpath=path.join(out_path,target_id)
    makedirs(outpath,exist_ok=True)

    
    

    
    clean_pdb = path.join(outpath,f"{target_id}_clean.pdb")
   

    structure = parsePDB(filepath)

    if not chain:
        chain = "A"
    
    
    protein   = structure.select(f"chain {chain} and not hetero")
    
    
    if protein is None:
        raise RuntimeError(f"Missing selected chain:{chain}")

    
    first_resnum = protein.getResnums()[0]
    
    if not altloc:
        altloc_groups = defaultdict(list)
        for atom in protein:
            alt = atom.getAltloc().strip()
            if alt:
                altloc_groups[alt].append(atom)


        
        if altloc_groups:
            altloc = max(
                altloc_groups.items(),
                key=lambda x: sum(a.getOccupancy() for a in x[1]) / len(x[1])
            )[0]
    
        
    writePDB(clean_pdb, protein)

    if repair:
        fixer = PDBFixer(filename=clean_pdb)
        fixer.findMissingResidues() #find and repalce missing residues in the file
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues() #find and replace non-standard resisdues
        if repair:
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()

        fix_pdb = path.join(outpath,f"{target_id}_fixed.pdb")
        with open(fix_pdb, 'w') as pdb_file:
            PDBFile.writeFile(fixer.topology, fixer.positions,pdb_file)
    else:
        fix_pdb = clean_pdb

    if forcefield:
        forcefield = forcefield.upper()
        amber_file = fix_pdb.replace("_fixed.pdb","_prot.pdb")
        pqr = fix_pdb.replace("_fixed.pdb","_prot.pqr")
        args = [
        f"--ff={forcefield}",   # forza campo
        f"--with-ph={pH}",
        "--ffout=AMBER",
        "--pdb-output",
        amber_file,          
        fix_pdb,
        pqr,
        ]

        out_file=fix_pdb.replace("_fixed.pdb",".pdb")
        pdb2pqr_main.run_pdb2pqr(args)


        mol = Chem.MolFromPDBFile(amber_file, removeHs=False)
        Chem.MolToPDBFile(mol, out_file)
    else:
        out_file = clean_pdb

    if cryst_ligand:
        cryst_ligand=cryst_ligand.upper()
        cryst_ligand_file = path.join(outpath,f"{target_id}_ligand.pdb")
        ligand = structure.select(f"resname {cryst_ligand} and not hydrogen and not water")
        if ligand:
            writePDB(cryst_ligand_file, ligand)
            cryst_mol = Chem.MolFromPDBFile(cryst_ligand_file,removeHs=False,sanitize=True)
            ligand_sdf_path = path.join(outpath,f"{target_id}_ligand.sdf")

            with Chem.SDWriter(ligand_sdf_path) as writer:
                writer.write(cryst_mol)

            return {"target_id":target_id,"file":out_file,"idx":first_resnum,"chain":chain,"altloc":altloc,"ligand":cryst_ligand_file,"ligand_sdf":ligand_sdf_path}
        else:
            raise RuntimeError(f"Passed crystallogrphic ligand residue name:{cryst_ligand} but is not found")
    else:
        return {"target_id":target_id,"file":out_file,"idx":first_resnum,"chain":chain,"altloc":altloc,"ligand":None}






class PipelineTemplate(ABC):
    def __init__(self,experiment_id:str,verbose:bool=True,run_name:str=None,source_dir:str="sources",ligand_dir:str="ligands",target_dir:str="targets",replicates:int=1,open_mp:bool=False):
        self.exp_id=experiment_id
        self.replicates=replicates
        base_path = "/home/screener/files"
        self.source_dir=path.join(base_path,source_dir)
        self.base_path=path.join(base_path,experiment_id) #workpath/expid/
        self.run_name=run_name
        if not run_name:
            self.work_path=path.join(self.base_path,f"run{randint(9999,99999)}")
        else:
            self.work_path=path.join(self.base_path,run_name)
        self.log=path.join(self.work_path,"screen.txt") #WIP
        self.maps_workpath=path.join(self.work_path,"maps")
        self.ligands_workpath=path.join(self.work_path,ligand_dir)
        self.targets_workpath=path.join(self.work_path,target_dir) #workpath/expid/ligands and targets
        self.result_workpath=path.join(self.work_path,"results")
        makedirs(self.source_dir,exist_ok=True)
        makedirs(path.join(self.source_dir,"target"),exist_ok=True)
        makedirs(self.ligands_workpath,exist_ok=True)
        makedirs(self.targets_workpath,exist_ok=True)
        makedirs(self.result_workpath,exist_ok=True)
        makedirs(self.maps_workpath,exist_ok=True)
        self.verbose=verbose
        self.docking_data=None
        self.target_data=None
        self.ligand_data=None #pandas log
        self.docking_report=None
        self.N_CONF_LIGAND=1
        try:
            mp.set_start_method("spawn", force=True)
        except RuntimeError:
            pass
        if open_mp:
            self.n_proc=int(getenv("SRUN_CPUS_PER_TASK"))
        else:
            self.n_proc=mp.cpu_count()-1

        self.docking_init = False  
    pass

    
    def load_ligands(self,filename:str,backend:str="rdkit",sanitize:bool=False):


        backend=backend.lower()
        #source_path=from_db
        source_path=path.join(self.source_dir,"ligands",filename)
        filetype=filename.split(".")[1]
        
        if backend == "rdkit":
            file_loader={"sdf":open_sdf,"sdf.gz":open_sdf,"csv":open_csv}
        else:
            raise RuntimeError(f"Unkown molecule backend: {backend}")
        
        loader=file_loader.get(filetype)

        if loader is None:
            supported=",".join(file_loader.keys())
            raise ValueError(f"Invalid file format:{filetype}\nSupported formats:{supported}")

        mols=loader(source_path,filetype,self.ligands_workpath)

        ligands= [{"id":m.GetProp("id"),"mol":m,"smiles":Chem.MolToSmiles(m)} for m in tqdm(mols,desc="Loading ligands")]
        
        with_duplicates = pd.DataFrame(ligands)

        no_dupes = with_duplicates.drop_duplicates(subset="smiles")
    
        if sanitize:
            tasks = ((row.smiles, row.id) for row in no_dupes.itertuples(index=False))
            
            n_workers = min(self.n_proc, len(tasks)) if len(tasks) else 1
            
            _fix_wrapper=partial(map_wrapper, func=fix_ligands)
            with mp.Pool(n_workers) as pool:
                results = []
                for r in tqdm(pool.imap_unordered(_fix_wrapper, tasks), desc="Sanitizing ligands"):
                    if r:
                        results.append(r)
            
            sanitized = pd.DataFrame(results) 
            print(f"Succefully loaded {sanitized.shape[0]},Discarded: {len(ligands)-sanitized.shape[0]}")
        else:
            sanitized = pd.DataFrame(ligands)

    

        self.__setattr__("ligand_data",sanitized) #join all results together
    pass
        
        
    def filter_ligands(self,filter_f:Callable,reject_indicator=None,*filter_args,**filter_kwargs):
        ligands = self.__getattr__("ligand_data")
        
        if ligands is None:
            raise RuntimeError("Missing ligand data")

        
        filter_f=partial(filter_f,*filter_args,**filter_kwargs)
        func_name=filter_f.func.__name__
        filtered_sets = []
        
        tasks = ((row.smiles, row.id) for row in ligands.itertuples(index=False))
        n_workers = min(self.n_proc, ligands.shape[0]) 
        #subdivide ligands list in chunks and elaborate them in mp
        
        with mp.Pool(n_workers) as pool:
            for r in tqdm(pool.imap_unordered(partial(map_wrapper, func=filter_f), tasks),
                total=ligands.shape[0],
                desc=f"Ligand filtering with {func_name}"):

                if not (r == reject_indicator):
                    filtered_sets.append(r)

       
        filtered=pd.DataFrame(filtered_sets)
        print(f"Passing: {filtered.shape[0]},Rejected: {ligands.shape[0]-filtered.shape[0]} ")   
        self.__setattr__("ligand_data",filtered) #join all results together

               
    pass

    
    def load_targets(self,source_target_files:Iterable=None,from_PDB:Iterable=None):
        

        print("Loading targets\n")
        if (source_target_files is None) and (from_PDB is None):
            raise Exception("Missing input target(s)")
        pass

        targets=[]

        if source_target_files:
            file_ids=[str((Path(file_path).stem)) for file_path in source_target_files]
            for i,_id in tqdm(enumerate(file_ids),desc="Loading target files"):
                targets.append({"target_id":_id,"file":source_target_files[i]})
        if from_PDB:
            dest_folder=path.join(self.source_dir,"target")
            pdb_entries = from_protein_data_bank(pdb_ids=from_PDB,destination_folder=dest_folder) #Download pdb files prom protetein databank and returns paths

        if pdb_entries:
            if not len(targets):
                targets= targets + pdb_entries
        
        self.__setattr__("target_data",pd.DataFrame(targets))

    
    def target_cleanup(self,keepWater:bool=False,prody_string:str=None,repair:bool=True,chunk_size:int=None,*args,**kwargs):
        targets = self.__getattr__("target_data")
        
        if not chunk_size:
            chunk_size=targets["target_id"].nunique()//self.n_proc
            if not chunk_size:
                chunk_size=targets["target_id"].nunique()

        clean_kwargs = {"out_path":self.targets_workpath,"keepWater":keepWater,"repair":repair}
        
        cleanup_f=partial(target_clean,**clean_kwargs)

        
        tasks = (
        (row.file, row.target_id, kwargs.get(f"{row.target_id}_altoc"), kwargs.get(f"{row.target_id}_chain"), kwargs.get(f"{row.target_id}_ligand"),kwargs.get("pH"),kwargs.get("forcefield"))
        for row in targets.itertuples(index=False)
        )

        cleaned_sets = []
        n_workers = min(self.n_proc, targets.shape[0]) 

        with mp.Pool(n_workers) as pool:
            for r in tqdm(pool.imap_unordered(partial(map_wrapper, func=cleanup_f), tasks),
                total=targets.shape[0],
                desc=f"Cleanup of {targets.shape[0]} target(s)"):
                
                cleaned_sets.append(r)
            
        
        cleaned=pd.DataFrame(cleaned_sets)
                
        

        self.__setattr__("target_data",cleaned) #join all results together
        print(f"\nSuccessfully cleaned {cleaned.shape[0]} target(s)\n")



    def check_poses(self,target_id:str,from_report:bool=False,redock:str=None,strict_steric_check:bool=False,threshold:float=2.0):

        lg = RDLogger.logger()
        lg.setLevel(RDLogger.ERROR)

        if from_report:
            report_file = self.__getattr__("docking_report")

            if not report_file:
                report_file=Path(f"{self.result_workpath}/{target_id}/{target_id}_{self.run_name}.csv")

                if not report_file.exists():
                    raise FileNotFoundError("Missing docking report")
            
            docking = pd.read_csv(report_file)
            docking["docked_sdf"] = docking["docked_sdf"].apply(ast.literal_eval)
            docking["energies"] = docking["energies"].apply(ast.literal_eval)
            
        else:
            docking = self.__getattr__("docking_data")

        busted_sets = []

        target_file = self.target_data.loc[target_id, "pdb_file"]
        bust_kwargs = {"target_id": target_id, "target_file": target_file, "strict_clash_check": strict_steric_check}
        buster = partial(bust_poses, **bust_kwargs)

        tasks = (
            (row.ligand_id, row.docked_sdf, row.energies, row.replicate)
            for row in docking.itertuples(index=False)
        )

        n_workers = min(self.n_proc, docking.shape[0]) if tasks else 1

        with mp.Pool(n_workers) as pool:
            for result in tqdm(pool.imap_unordered(partial(map_wrapper,func=buster), tasks),total=docking.shape[0],desc=f"Checking docking poses for {target_id}"):
                busted_sets.append(result)


        busted_mols = pd.DataFrame(busted_sets)

        busted_mols.to_csv(f"{self.result_workpath}/{target_id}/{target_id}_{self.run_name}_busted.csv",encoding="utf-8",sep=",",index=False)

        replicated_sets=[]

        busted_group=busted_mols.groupby("ligand_id").agg({
            "file": tuple,
            "energy": tuple
            }).reset_index()

        
        replicated_sets = []

        repro_checker = partial(check_producibility, target_id=target_id, threshold=threshold)

        tasks = (
            (row.ligand_id, row.file, row.energy)
            for row in busted_group.itertuples(index=False)
        )

        n_workers = min(self.n_proc, busted_group.shape[0]) if tasks else 1

        with mp.Pool(n_workers) as pool:
            for result in tqdm(pool.imap_unordered(partial(map_wrapper,func=repro_checker), tasks),total=busted_group.shape[0],desc=f"Checking reproducibility of poses for {target_id}"):
                replicated_sets.append(result)
            
        replicated = pd.DataFrame(replicated_sets)
        replicated.to_csv(f"{self.result_workpath}/{target_id}/{target_id}_{self.run_name}_reproduced.csv",encoding="utf-8",sep=",",index=False)
        
        if redock:
            replicated.set_index("ligand_id",inplace=True)
            print(f"\nTesting redocking on {redock} for target {target_id}...")

            redock_status = "No"
            
            cryst_mol = SDMolSupplier(self.target_data.loc[target_id,"ligand_sdf"])[0]
            docked_mol = SDMolSupplier(replicated.loc[redock,"file"])[0]

            SanitizeMol(cryst_mol)
            SanitizeMol(docked_mol)

            o3a = AllChem.GetO3A(docked_mol, cryst_mol)
            rmsd = o3a.Align()

            if rmsd <= threshold:
                redock_status="Yes"

            print(redock_status)

            for ligand in replicated_sets:
                if ligand["ligand_id"] != redock:
                    ligand["redocked"] = "Not in structure"
                else:
                    ligand["redocked"] = redock_status

            
            report_df=pd.DataFrame(replicated_sets)

            self.__setattr__("report",report_df)
            
            report_df.to_csv(f"{self.result_workpath}/{target_id}/{target_id}_{self.run_name}_redock.csv",encoding="utf-8",sep=",",index=False)
        else:
            self.__setattr__("report",replicated)


    def calc_metrics(self,parameters,target_id:str,temp:float=37,energy_unit ="kcal",round_digit:int=5):
        
        energy_unit=energy_unit.lower()
        supported_energy_unit=("kcal","kj")

        if energy_unit not in supported_energy_unit:
            raise ValueError("Unsupported energy unit:{}\nSupported:{}".format(energy_unit,supported_energy_unit.join(",")))

        metrics = self.__getattr__("report")
        metrics.reset_index(drop=False)
        
 
        


        metrics_sets = []

        tasks = (
            (row.ligand_id,
            row.file,
            target_id,
            row.energy,
            temp,
            energy_unit,
            round_digit,
            parameters)
            for row in metrics.itertuples(index=False)
        )

        n_workers = min(self.n_proc, metrics.shape[0]) if tasks else 1

        with mp.Pool(n_workers) as pool:
            for result in tqdm(pool.imap_unordered(partial(map_wrapper,func=calc_parameters), tasks),total=metrics.shape[0],desc=f"Calculating parameters for {target_id}"):
                metrics_sets.append(result)

        metrics_df = pd.DataFrame(metrics_sets)

        full_report = pd.merge(metrics_df,metrics,on=["ligand_id","target_id"],how="left")
        self.__setattr__("full_report",full_report)
        full_report.to_csv(f"{self.result_workpath}/{target_id}/{target_id}_{self.run_name}_final.csv",encoding="utf-8",sep=",",index=False)
       

        


            


    @abstractmethod
    def init_docking(self,*args,**kwargs):
        ...

    @abstractmethod
    def optimize_ligand(self,*args,**kwargs):
        ...

    @abstractmethod
    def process_targets(self,*args,**kwargs):
        ...
    
    @abstractmethod
    def process_ligands(self,*args,**kwargs):
        ...

    
    






    """
    def screen_single_target(self,target_id:str,report_fn:str=None,dock_kwargs:dict=None):

        if not report_fn:
            report_fn = f"{target_id}_{self.run_name}"





        ligand_ids=self.ligand_data.index
        replicates={}
        for j in range(1,self.replicates+1):
            agent.reset_seed()
            print(f"Replicate:{j}")
            ligand_energies={"target id":target_id}

            avg_time = {}
            for ligand in ligand_ids:
                times = []
                fn = coupled_fn[ligand]
                try:
                    #ligand_mol=self.ligand_data.loc[ligand]["mol"]
                    n_conf=self.N_CONF_LIGAND
                    result_path=f"{self.result_workpath}/{ligand}"
                    makedirs(result_path,exist_ok=True)
                    start = time.perf_counter()
                    if dock_kwargs:
                        base_energy=agent.__dock__(ligand_pdbqt=self.ligand_data.loc[ligand]["base"],ligand_id=ligand,target_id=target_id,coupled_fn=fn,save_path=result_path,replicate_id=j,target_pdbqt=target,**dock_kwargs)
                    else:
                        base_energy=agent.__dock__(ligand_pdbqt=self.ligand_data.loc[ligand]["base"],ligand_id=ligand,target_id=target_id,coupled_fn=fn,save_path=result_path,replicate_id=j,target_pdbqt=target)
                    end = time.perf_counter()

                    times.append(end-start)
                    ligand_energies[ligand]=base_energy
                    if n_conf != -1:
                        conf_energy=[]
                        #ligand_energies[f"base_{ligand}"]=base_energy
                        confs_pdbqt=self.ligand_data.loc[ligand]["confs"]
                        for i in range(n_conf):
                            start = time.perf_counter()
                            if dock_kwargs:
                                energy=agent.__dock__(ligand_pdbqt=confs_pdbqt[i],ligand_id=ligand,target_id=target_id,coupled_fn=fn,conf_id=i,save_path=result_path,replicate_id=j,target_pdbqt=target,**dock_kwargs)
                            else:
                                energy=agent.__dock__(ligand_pdbqt=confs_pdbqt[i],ligand_id=ligand,target_id=target_id,coupled_fn=fn,conf_id=i,save_path=result_path,replicate_id=j,target_pdbqt=target)
                            conf_energy.append(energy)
                            end = time.perf_counter()
                            times.append(end-start)
                            ligand_energies[f"conf_{i}_{ligand}"]=energy
                        pass
                    pass
                    avg_time[f"total time {ligand}"] = sum(times)
                except Exception as e:
                    print(ligand)
                    print(e)
                    continue
            pass
            #print(ligand_energies)
            replicates[f"replicate {j}"]=ligand_energies | avg_time
            self.docking_data=pd.DataFrame.from_dict(replicates,orient="index")
            #print(self.docking_data)
            self.docking_data.to_csv(f"{self.result_workpath}/{report_fn}",encoding="utf-8",sep=",")
            #print(f"Time Elapsed{start-end},time per opearation avg {(start-end)/(n_conf*j*i)}")
        pass
    pass


    screen_multi_target: Screens multiple target proteins using the same docking agent
    Paramenters:
        map_modes: dict containing map modes for each target e.g. for proteins aaa and bbb the dicts should be {"aaa":"o4a","bbb":"flex"}
        agent: Callable docking agent class
        multinode: toggle for multi node computing, each node works on a single target
    """


    def screeen_multi_target(map_modes:dict,agent:Callable,multinode:bool=False):
        targets=self.target_data["target_id"]
        if not multinode:
            for target in targets:
                Screener.screen_single_target(map_mode=map_modes[target],docking_agent=agent)
        else:
            raise NotImplementedError
    


                

    

    def __getattr__(self, name: str):
        return self.__dict__[f"_{name}"]

    def __setattr__(self, name:str, value):
        self.__dict__[f"_{name}"] = value

    @property
    def ligands(self):
        return self.ligand_data.copy(deep=True)

    @property
    def targets(self):
        return self.target_data.copy(deep=True)

    @property
    def docking(self):
        return self.docking_data.copy(deep=True)