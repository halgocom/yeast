from abc import ABC,abstractmethod
from dataclasses import dataclass
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from prody import parsePDB,saveAtoms,writePDB,AtomGroup
from os import getcwd,listdir,path,mkdir,makedirs,getenv
from typing import Iterable
from loader import read_input
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
from rdkit import RDLogger,Chem
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




def bust_poses(ligand_id:str,
dock_results:list,
dock_energies:list,
replicate:int,
target_id:str,
target_file:str,
box_center:Iterable,
box_size:Iterable,
strict_clash_check:bool=False,
box_threshold:float=0,
channel_pocket:dict=None,
score_threshold:float=-5.0):


    
    box_lower = box_center - box_size / 2
    box_upper = box_center + box_size / 2
    fractions = np.empty(len(dock_results), dtype=float)

    if channel_pocket:
        pocket_center = np.array(channel_pocket[0])
        pocket_size = np.array(channel_pocket[1])

        pocket_lower = pocket_center - pocket_size / 2
        pocket_upper = pocket_center + pocket_size / 2




    for i,pose in enumerate(dock_results):
        mol = SDMolSupplier(pose)[0]
        pos = mol.GetConformer().GetPositions()
        
        #find heavy atom indicies
        heavy = np.array(
        [atom.GetAtomicNum() > 1 for atom in mol.GetAtoms()],
        dtype=bool
        )

        #keep only heavy atom
        pos = pos[heavy]

        in_box = (pos >= box_lower) & (pos <= box_upper)
        if channel_pocket:
            in_pocket = (pos >= pocket_lower) & (pos <= pocket_upper)
            pass_condition = in_box & in_pocket
        else:
            pass_condition = in_box

            
        fractions[i] = np.all(pass_condition,axis=1).mean()

    dock_energies = np.asarray(dock_energies)
    mask = (fractions >= box_threshold) & (dock_energies <= score_threshold)

    dock_results = np.asarray(dock_results)[mask]
    dock_energies = np.asarray(dock_energies)[mask]
    fractions = fractions[mask]



    if len(dock_results) == 0 or np.all(dock_energies > score_threshold):
        return None


        

    
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

    if best_confs.empty:
        return None
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


    




def target_clean(filepath:str,
target_id:str,
altloc:str=None,
chain:str=None,
cryst_ligand:str=None,
pH:float=7.4,
forcefield:str=None,
keepWater:bool=False,
repair:bool=True,
out_path:str=None):
    
    
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

            lig_mol = Chem.SDMolSupplier(ligand_sdf_path)[0]
            center,size=get_box(lig_mol)

            with open(path.join(outpath,f"{target_id}_pocket.pdb"),"w") as f:
                f.write(box_to_pdb_string(center,size))
            return {"target_id":target_id,"file":out_file,"idx":first_resnum,"chain":chain,"altloc":altloc,"ligand":cryst_ligand_file,"ligand_sdf":ligand_sdf_path,"ligand_box":(center,size)}
        else:
            raise RuntimeError(f"Passed crystallogrphic ligand residue name:{cryst_ligand} but is not found")
    else:
        return {"target_id":target_id,"file":out_file,"idx":first_resnum,"chain":chain,"altloc":altloc,"ligand":None,"ligand_box":None}






class PipelineTemplate(ABC):
    def __init__(self,experiment_id:str,verbose:bool=True,run_name:str=None,source_dir:str="sources",ligand_dir:str="ligands",target_dir:str="targets",replicates:int=1,open_mp:bool=False):
        self.exp_id=experiment_id
        self.replicates=replicates
        base_path = "/home/screener/files"
        self.source_dir=path.join(base_path,source_dir)
        self.base_path=path.join(base_path,experiment_id) #workpath/expid/
        
        if not run_name:
            self.run_name=f"run{randint(9999,99999)}"
            self.work_path=path.join(self.base_path,self.run_name)
        else:
            self.run_name=run_name
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

    def load_ligands(self,filename:str, sanitize:bool=True, sdf_confs:bool=False, sep:str='\t'):

        source_path=path.join(self.source_dir,"ligands",filename)

        suppl = read_input(fname=source_path,id_field_name="ligand_id",sanitize=sanitize,sdf_confs=sdf_confs)

        ligands = [{"id":_id,"mol":mol,"smiles":smiles} for mol,_id,smiles in tqdm(suppl,desc="Loading ligands")]

        
        ligands_df = pd.DataFrame(ligands)

        self.__setattr__("ligand_data",ligands_df) 
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

    
    #def load_targets(self,source_target_files:Iterable=None,from_PDB:Iterable=None):
    def load_targets(self,*target_names):
        
        if len(target_names) ==0:
            raise Exception("Missing input target(s)")
        

        src_folder=Path(self.source_dir)/"target"
        targets=[]

        for target in tqdm(target_names,desc="Loading targets",total = len(target_names),position=0):

            if not isinstance(target,str):#filter only targets with string names
                continue

            file_path = src_folder / target

            if file_path.is_file(): #if file exits loads it directly
                filename = target.split(".")[0]
                targets.append({"target_id":filename,"file":str(file_path)})
            else: #tries to download it
                try:
                    download = from_protein_data_bank(pdb_id=target,destination_folder=src_folder)
                except Exception as e:
                    print(f"Target {target} errored while downloading the Protein Data Bank with error:\n{e}")
                    continue

                targets.append(download)
        pass

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
        (row.file, row.target_id, 
        kwargs[row.target_id].get("altoc"), 
        kwargs[row.target_id].get("chain"),
        kwargs[row.target_id].get("ligand"),
        kwargs[row.target_id].get("pH"),
        kwargs[row.target_id].get("forcefield"))
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



    def check_poses(self,target_id:str,from_report:bool=False,redock:str=None,strict_steric_check:bool=False,box_threshold:float=0,rsmd_threshold:float=2.0,channel_pocket:bool=False):

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
        

        box_center,box_size = self.target_data.loc[target_id,["box_center","box_size"]]
        box_center,box_size = np.array(box_center,dtype=float),np.array(box_size,dtype=float)
        if channel_pocket:
            pocket= self.target_data.loc[target_id,"ligand_box"]
            bust_kwargs = {"target_id": target_id, 
            "target_file": target_file, 
            "strict_clash_check": strict_steric_check,
            "box_center":box_center,
            "box_size":box_size,
            "box_threshold":box_threshold,
            "channel_pocket":pocket}
        else:
            bust_kwargs = {"target_id": target_id, 
            "target_file": target_file, 
            "strict_clash_check": strict_steric_check,
            "box_center":box_center,
            "box_size":box_size,
            "box_threshold":box_threshold}

            
        buster = partial(bust_poses, **bust_kwargs)

        tasks = (
            (row.ligand_id, row.docked_sdf, row.energies, row.replicate)
            for row in docking.itertuples(index=False)
        )

        n_workers = min(self.n_proc, docking.shape[0]) if tasks else 1

        with mp.Pool(n_workers) as pool:
            for result in tqdm(pool.imap_unordered(partial(map_wrapper,func=buster), tasks),total=docking.shape[0],desc=f"Checking docking poses for {target_id}"):
                if result:
                    busted_sets.append(result)


        busted_mols = pd.DataFrame(busted_sets)

        busted_mols.to_csv(f"{self.result_workpath}/{target_id}/{target_id}_{self.run_name}_busted.csv",encoding="utf-8",sep=",",index=False)

        replicated_sets=[]

        busted_group=busted_mols.groupby("ligand_id").agg({
            "file": tuple,
            "energy": tuple
            }).reset_index()

        
        replicated_sets = []

        repro_checker = partial(check_producibility, target_id=target_id, threshold=rsmd_threshold)

        
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

            try:
                docked_mol = SDMolSupplier(replicated.loc[redock,"file"])[0]
            except Exception as e:
                print(f"Missing docked crystallografic molecule {redock}")
                self.__setattr__("report",replicated)
                return None

            SanitizeMol(cryst_mol)
            SanitizeMol(docked_mol)

            o3a = AllChem.GetO3A(docked_mol, cryst_mol)
            rmsd = o3a.Align()

            if rmsd <= rsmd_threshold:
                redock_status="Yes"

            print(redock_status)

            for ligand in replicated_sets:
                if ligand["ligand_id"] != redock:
                    ligand["redocked"] = "Not in structure"
                else:
                    ligand["redocked"] = redock_status

            
            report_df=pd.DataFrame(replicated_sets)

            report_df["norm_score"] = np.round(report_df["energy"] / (report_df["energy"].min()/100),2)
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
        self.__setattr__("report",full_report)
        full_report.to_csv(f"{self.result_workpath}/{target_id}/{target_id}_{self.run_name}_final.csv",encoding="utf-8",sep=",",index=False)
       

        


            


    @abstractmethod
    def init_docking(self,*args,**kwargs):
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