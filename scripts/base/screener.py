from AutoDockScreen import VinaScreen
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from random import randint
from meeko import PDBQTReceptor
from meeko import PDBQTMolecule
from meeko import ResidueChemTemplates
from meeko import Polymer
from prody import parsePDB,checkNonstandardResidues, saveAtoms,writePDB,calcCenter
from meeko.gridbox import calc_box
from rdkit.Geometry import Point3D
from os import getcwd,listdir,path,mkdir,makedirs,getenv
from typing import Iterable
import pandas as pd
import numpy as np
from pathlib import Path
from collections.abc import Callable
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import multiprocessing as mp
from multiprocessing import Pool
import csv
import subprocess
import sys
import time
import gzip
import random
import string 
import shutil
from molscrub import Scrub
from urllib.request import urlretrieve
from pdbfixer import PDBFixer
from openmm.app import PDBFile
#from AutoDockTools import MoleculePreparation as molprep
from reduce2_interface import Reduce2App
from dotenv import load_dotenv 
from lib import *
#from mmtbx.programs import reduce2
#from iotbx.cli_parser import run_program

#import AutoDockTools as ADT
"""
Struttura Pipeline:
Inserimento
    Target
        da db
        Locale *
        internet 
    Ligando
        Da db
        Locale *
        Internet
    Impostazioni
        json
        [Creare un file statico tipo docker per farlo in automatico]

Preparazione
    Target
        Import
        per ogni target
            crea la sc
    ligando
        Import
        per ogni ligando
            Ottenimento numero configurazioni nconf
            Genenrazione conf
            salvetaggio di ogni conf in un file pdbqt (./lig1/conf1.pdbqt)
            genenrare il gridbox per ogni conformazione
            dopo aver salvato ogni configurazione se serve idratare 


Generazione file di campo
    Autogrid
        da creare ad ogni ligando -> lunghissimo si potrebbe threddare -> lo facciamo prima 
        Gestione configurazioni multiple -> Ogni configurazione dello stesso ligando ha le stesse mappe
    usando autoligand si può trovare le migliori dimensioni della box e centro e per ricoprire [da fare]
        
Docking
    Vina
    autodock 4
    Protenix

Analisi
    

"""

"""
TODO:
[X] Creare lista/tupla per percorsi ligandi per successivo screening
[X] Creare metodo preparazione target
Usare allign per risolvere il prblema del ligando cristallo grafico estrai converti in pdb e apri come rdkit mol
Creare metodo preparazione target gennato
Creare metodo istanza docking (therad o multiproces mo vediamo) single target (multitraget per ora improbabile)
Metodo di logging  e Data collection  con confornto
"""
"""
Librerie:
    Meeko: Interfraccia per autodock VINA
    rdkit
"""


"""
    Constructor parameters

    ligand_dir:str
        default:"ligands",folder where ligand(s) files are saved
    target_dir:str
        default:"targets",folder where target protein(s) files are saved
    verbose:bool
        default:True enables verbose output
    work_path:str
        Absolute path where all docking ligand and target source files are located
    verbose:bool
        Enables verbose output
    
    Class Attributes

    work_path: Same as parameter
    ligands_workpath: Absolute path to ligand folder
    targets_workpath: Absolute path to target folder 
    verbose: same as paramenter
    ligands: ligand(s) name(s) cache
    ligands_pdbqt: ligand(s) pdbqt cache
    data :pandas dataframe for data collection #si potrebbe creare un oggetto apposta da ereditare se necessario
"""





class Screener:
    def __init__(self,experiment_id:str,verbose:bool=True,run_name:str=None,source_dir:str="sources",ligand_dir:str="ligands",target_dir:str="targets",replicates:int=1):
        self.exp_id=experiment_id
        self.replicates=replicates
        load_dotenv(dotenv_path="/home/screener/env/.env")
        base_path=getenv("FILE_PATH")
        self.source_dir=path.join(base_path,source_dir)
        self.base_path=path.join(base_path,experiment_id) #workpath/expid/
        if not run_name:
            self.work_path=path.join(self.base_path,f"run{randint(9999,99999)}")
        else:
            self.work_path=path.join(self.base_path,run_name)
        self.log=path.join(self.work_path,"screen.txt") #WIP
        self.maps_workpath=path.join(self.work_path,"maps")
        self.ligands_workpath=path.join(self.base_path,ligand_dir)
        self.targets_workpath=path.join(self.base_path,target_dir) #workpath/expid/ligands and targets
        self.result_workpath=path.join(self.work_path,"results")
        #makedirs(path.join(self.source_dir,"ligands"),exist_ok=True)
        makedirs(self.source_dir,exist_ok=True)
        makedirs(path.join(self.source_dir,"target"),exist_ok=True)
        makedirs(self.ligands_workpath,exist_ok=True)
        makedirs(self.targets_workpath,exist_ok=True)
        makedirs(self.result_workpath,exist_ok=True)
        makedirs(self.maps_workpath,exist_ok=True)
        #print(f"Directory created:{self.ligands_workpath}\n{self.targets_workpath}")
        self.verbose=verbose
        #self.ligands_pdbqt={}
        #self.ligands=[] #Ligands()
        self.docking_data=None
        self.target_data=None
        self.ligand_data=None #pandas log
        self.N_CONF_LIGAND=0 
        self.bin_dir=path.dirname(path.abspath(__file__)) #binary folder path, for _script
        self.python_path=str(sys.executable) #python executable path, for _script
    pass

    @classmethod
    def from_config(cls,config,json=False):

        pass

    #WIP
   


    """
    prepare_ligand
        ligand_id:str process id
        mol:RDKit mol object
        n_conf:int  Number of conformation to generate
        hydrate:bool default:False activate hydrated protocol
        prepation_config:dict default:None meeko preparator configuration dict
        output:str default: file Specify needed output if pdbqt needed use file if path is needed use path
    """

    def mk_prepare_ligand(self,mol,receptor_id:str,ligand_id:str,n_conf:int,box_mode:str="geometric",hydrate:bool=False,output:str="file",preparation_config:dict=None,etkdg_seed:int=65):

        
        print(f"Process {ligand_id} started")
        mol=Chem.AddHs(mol)
        ligand_dir=path.join(self.ligands_workpath,ligand_id)
        
        # Molecule preparation from Meeko
        if preparation_config is not None:
            preparator = MoleculePreparation.from_config(preparation_config) #load config
        else:   
            preparator = MoleculePreparation(hydrate=hydrate) #hydrate protocol

        makedirs(ligand_dir,exist_ok=True)
        
        
        
        center,size=self.get_box(mol=mol, box_mode=box_mode)
        center=np.array(center)
        size=np.array(size)



        ligand=preparator.prepare(mol)
        ligand_pdbqt = PDBQTWriterLegacy.write_string(ligand[0])
        ligand_pdbqt=self.write_pdbqt(dir=ligand_dir,filename=ligand_id,pdbqt=ligand_pdbqt) #save base ligand pdbqt
        if output == "path":
            ligand_pdbqt=f"{ path.join(ligand_dir, f"{ligand_id}.pdbqt") }"
        
        # Generate conformations with 3D atom positions
        if n_conf >1:
            confs=[]
            sizes=[]
            
            #makedirs(conf_dir, exist_ok=True)

            params = AllChem.ETKDGv3() #experimental torsion knowledge distance geometry (ETKDG) method lo usiamo perchè quello di ottimizzazione stanstard di rdkit non è compatbile con mp
            params.useSmallRingTorsions = True
            params.useMacrocycleTorsions =True
            params.randomSeed = etkdg_seed
            params.trackFailures = True
            params.enforceChirality = True
            params.clearConfs = True
            conf_ids = AllChem.EmbedMultipleConfs(mol,numConfs=n_conf,params=params)
            rdMolAlign.AlignMolConformers(mol)
            for conf_id in conf_ids:
                
                Screener.coords_align(mol,center,conf_id)
                _,conf_size=self.get_box(mol=mol,box_mode=box_mode,conf_id=conf_id)
                sizes.append(conf_size)
                #centers.append(center)
                
                if hydrate: #addictional flags for hydration      
                    preparator = MoleculePreparation(hydrate=hydrate) #hydrate protocol
                conf=preparator.prepare(mol,conformer_id=conf_id)
                
                conf_pdbqt = PDBQTWriterLegacy.write_string(conf[0])
                self.write_pdbqt(dir=ligand_dir,filename=f"{ligand_id}_{conf_id}",pdbqt=conf_pdbqt) #save conf pdbqt
                
                if output == "file":
                    confs.append(conf_pdbqt[0])
                else:
                    confs.append(f"{ path.join(ligand_dir, f"{ligand_id}_{conf_id}.pdbqt") }")
            sizes=np.array(sizes)
            opt_size=sizes.max(axis=0)
        elif n_conf == 1:
            return {"ligand_id":f"{ligand_id}","mol":mol,"base":ligand_pdbqt,"box_center":center,"box_size":size}
        else:
            raise ValueError(f"Invalid number of configurations :{n_conf}")
        
        return {"ligand_id":f"{ligand_id}","mol":mol,"base":ligand_pdbqt,"confs":tuple(confs),"box_center":center,"box_size":opt_size}
    pass

    """
    process_ligand
        filename:str
        n_conf:int default:1 Number of conformation to genenrate
        hydrate:bool default:False activate hydrated protocol
        prepation_config:dict default:None meeko prepaprator configuration dict
    """
    
    """
    6.4–6.8 sano Metabolic changes during prostate cancer development and progression
    6.1-6.5 malato The pH of prostatic fluid: a reappraisal and therapeutic implications.
    """
    def process_ligands(self,
        filename:str,
        n_conf:int=5,
        n_proc:int=None,
        ph:dict={"low":7.4,"high":7.4},
        preparation:Callable=mk_prepare_ligand,
        tasks:list=None,
        box_mode:str="geometric",
        hydrate=False,
        addH:bool=True,
        output:str="file"):
        
        print("Start ligand preparation")
        start = time.perf_counter()
        self.N_CONF_LIGAND=n_conf
        if n_proc is None:
            n_proc=mp.cpu_count()-2 #use max cores #remove -2 for HPC
        
    
        #source_path=from_db
        source_path=path.join(self.source_dir,"ligands",filename)
        filetype=filename.split(".")[1]
        file_loader={"sdf":Screener.open_sdf,"sdf.gz":Screener.open_sdf,"csv":Screener.open_csv}
        
        loader=file_loader.get(filetype)

        if loader is None:
            supported=",".join(file_loader.keys())
            raise ValueError(f"Invalid file format:{filetype}\nSupported formats:{supported}")

        mols=loader(source_path,filetype,self.ligands_workpath)
        
        scrub=Scrub(ph_low=ph["low"],ph_high=ph["high"]) #Protonation agent
        #scrubbed_mols=tuple([scrub(mol)[0] for mol in mols]) #protoned molecules

        if addH:
            scrubbed_mols=[]
            for mol in mols:
                print(mol.GetProp("id"))
                scrubbed_mols.append(scrub(mol)[0])

            scrubbed_mols=tuple(scrubbed_mols)
        else:
            scrubbed_mols=tuple(mols)

        for receptor_id in self.target_data.index:
            box_origin=self.target_data.loc[receptor_id,"box_type"]
            if box_origin == "LIG":
                box_mol,receptor_center,receptor_size = self.target_data.loc[receptor_id,["box_mol","box_center","box_size"]].to_numpy()
                receptor_box_center,receptor_box_size = np.array(receptor_center), np.array(receptor_size)
                
                box_mol_conf=box_mol.GetConformer()
                
                for mol in scrubbed_mols:
                    mol_conf=mol.GetConformer()
                    
                    try: #try align ligand to cristallographic ligand
                        rdMolAlign.AlignMol(prbMol=mol_conf,refMol=box_mol_conf,maxIters=200) 
                    except Exception:
                        pass
            else:
                receptor_box_center,receptor_box_size = self.target_data.loc[receptor_id,["box_center","box_size"]].to_numpy()

            aligned_mols=tuple([Screener.coords_align(mol,receptor_box_center) for mol in scrubbed_mols])

            
            #ids=tuple([mol.GetProp("id") for mol in aligned_mols])


            ligands={}
            for mol  in aligned_mols:
                lid = mol.GetProp("id")
                ligands[lid] = mol


            if tasks is  None:
                tasks = [(self,ligands[ligand_id],receptor_id,ligand_id, n_conf, box_mode, hydrate, output) for ligand_id in ligands.keys()]  #calculate task paramenters
            
            with Pool(n_proc) as pool:
                results=pool.starmap(preparation, tasks)

            #sostituire questo con dataframe in shared memory
            data={}
            for column in results[0].keys():
                data[column] = tuple(result[column] for result in results)

            self.ligand_data=pd.DataFrame(data)
            self.ligand_data=self.ligand_data.set_index("ligand_id")

      
        end = time.perf_counter()
        
        print(f"Terminated in {end-start:.3f}s,avg:{(end-start)/len(results):.3f} s/ligand")
    pass

    """
        prepare_target: prepares single pdbqt receptor file
            file_path:string, pdb/pdbqt/mmCIF absolute filepath
            receptor_id:string, receptor_id extracted from filename or use PDB id
            hydrate: bool, hydrate flag
            flexres: list,list of flexible residues string formatted as "chain_name:residue_name_residue_position 
            e.g. to target alanine residue in the 512th position on chain A use 'A:ALA:512' "
            backend: backend to use, defaults to prody
    """


    def mk_prepare_target(self,
        box:Iterable,
        file_path:str,
        receptor_id:str,
        hydrate:bool=False,
        charge_type="gasteiger",
        flexres:list=None,
        config:dict=None,
        backend:str="prody",
        ph:float=7.0
        ):


        if config is not None:
            preparation_config=config.get("preparation_config")
            blunt_ends=config.get("blunt_ends")
            wanted_altloc=config.get("wanted_altloc")
            delete_residues=config.get("delete_residues")
            allow_bad_res=config.get("allow_bad_res")
            set_template=config.get("set_template")
            charge_type=config.get("charge_type")
        file_format=file_path.split(".")[-1]


        #fare in modo che reduce non dia più errore trasferendo questo il pezzo della scelta del converter  in process

        
        base_pdb_path = path.join(self.source_dir,"target",f"{receptor_id}.pdb")
        receptor_pdb_path = path.join(self.targets_workpath,receptor_id,f"{receptor_id}_H.pdb")
        outpath=path.join(self.targets_workpath,receptor_id,f"{receptor_id}.pdbqt")
        
        
        if isinstance(box,str):
            box_type,origin=box.split(":")
            auto_box={"LIG":Screener.cligand_box,"RES":Screener.residue_box}
            box_calculator=auto_box.get(box_type)
            if box_calculator is None:
                raise ValueError(f"Invalid box aurgument:{box_type}")
            box=box_calculator(self,base_pdb_path,origin)
        else:
            box={"box_type":"USER","box_center":box[0],"box_size":box[1]}
        
        
        
 
            
        
        #set target preparation parameters
        templates = ResidueChemTemplates.create_from_defaults() #create from defaults for now
        if config is not None:
            mk_prep = MoleculePreparation.from_config(config)
        else:
            mk_prep = MoleculePreparation(hydrate=hydrate,charge_model=charge_type)
    
        #load pdb with H
        if backend.lower() == "prody":
            target=self.prody_load(receptor_pdb_path)

            polymer=Polymer.from_prody(
            target,
            templates,  # residue_templates, padders, ambiguous,
            mk_prep,
            #set_template,
            #delete_residues,
            allow_bad_res=True,
            #blunt_ends=True,
            #wanted_altloc=wanted_altloc,
            default_altloc="A"
        )
    
        elif backend.lower() == "file":
            with open(receptor_pdb_path,"r") as pdb_file:
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
            default_altloc="A",
        )
        else:
            raise ValueError(f"Invalid back end:{backend}")


        #pdb_string=mmCIF_to_pdb()
        
        
        
        pdbqt_tuple = PDBQTWriterLegacy.write_from_polymer(polymer)
        pdbqt_path=f"{self.targets_workpath}/{receptor_id}/{receptor_id}.pdbqt"
        pdbqt=Screener.write_pdbqt(dir=f"{self.targets_workpath}/{receptor_id}",filename=f"{receptor_id}",pdbqt=pdbqt_tuple)
        receptor=PDBQTReceptor.from_pdbqt_filename(outpath)

        #print(Screener.__get__("cryst_lig_center"))
        return {"id":receptor_id,"path":pdbqt_path} | box #return dict with receptor and box specs
        

        if flexres is not None:
            flex_residues = set()
            for string in flexres:
                for res_id in parse_residue_string(string):
                    flex_residues.add(res_id)
            for res_id in flex_residues:
                receptor.flexibilize_sidechain(res_id, mk_prep)
        
        

            pdbqt_tuple = PDBQTWriterLegacy.write_from_polymer(polymer)
            rigid_pdbqt, flex_pdbqt_dict = pdbqt_tuple
            rigid=Screener.write_pdbqt(dir=f"{self.targets_workpath}/{receptor_id}",filename=f"{receptor_id}_rigid",pdbqt=rigid_pdbqt)
            if flexres:
                flex=Screener.write_pdbqt(dir=f"{self.targets_workpath}/{receptor_id}",filename=f"{receptor_id}_flex",pdbqt=flex_pdbqt)
            pass
            return {"id":receptor_id,"pdbqt":rigid_pdbqt,"flex":flex_pdbqt_dict,"flexres":flexres} | box #return dict with receptor and box specs


        """
            Per risolver il proplema del multiprocessing con pandas esisite un sisitema anche per il caricemnto dei plugin,per i divertsi tipi di box (quella base e quella ottimizzata del paper) 
            e il numero metti una variaibile "unapertutti(questa autimatica e corrispondente al recettore per flexres)/unaperognuno"
        """
        
    pass



    """
    process_targets: Parallelized processing of target proteins
        source_receptor_files:list of strings containing target pdb filenames in source folder
        from_PDB: list of string containing Protein DataBase id of target(s)
        hydrate:bool toggle hydrate protocol
        preparation: preparation function defaults to one given in its class
        n_proc:int number of processors used in mp if None default to max processors leaving 2 for OS
        flexres: dict of  lists of strings containing name and position of flexible residue for each taget  eg "A:ALA:512" targets the alanine residue in the 512 place in chain A 
        map_mode:str Specify if need to calculate one map for each ligand or use one map for all
        -Optional arguments
            prody_string: dict formatted as {"target_id":"prody string"} if using prody as pdb backend this defaults to only select the protein backbone structure and water (if specified using hydrate)
    
    
    
    
    """


    def process_targets(self,
        box:Iterable,
        source_receptor_files:Iterable=None,
        from_PDB:Iterable=None,
        n_proc:int=None,
        hydrate:bool=False,
        preparation:Callable=mk_prepare_target,
        flexres:dict=None,
        map_mode:str="o4e",
        charge_type="geisteger",
        config:dict=None,
        ph:float=7.0,
        backend:str="prody",
        *args,
        **kwargs
        ):

        start = time.perf_counter()
        if (source_receptor_files is None) and (from_PDB is None):
            raise Exception("Missing input receptor(s)")
        pass

        if n_proc is None:
            n_proc=mp.cpu_count()-2 #use max cores #remove -2 for HPC
        
        
        receptors={}

        if source_receptor_files:
            file_ids=[str((Path(file_path).stem)) for file_path in source_receptor_files]
            for i,_id in enumerate(file_ids):
                receptors[_id]=source_receptor_files[i]
        if from_PDB:
            pdb_entries = Screener.from_protein_data_bank(self,from_PDB) #Download pdb files prom protetein databank and returns paths

        if pdb_entries:
            receptors= receptors | pdb_entries

                

        #clean up and add hydrogens to receptor before preparation


        pdb_cleaner={"pdbfix":Screener.pdb_fix,"prody":Screener.prody_fix}
        PDB=pdb_cleaner.get(backend)

        if PDB is None:
            raise ValueError(f"Invalid backend\nSupported: pdbfix,prody")
        #box={}
        
        tasks=[]
        for _id in receptors.keys():
            makedirs(path.join(self.targets_workpath,_id),exist_ok=True)
            receptors[_id] = PDB(receptors[_id],self.targets_workpath,_id,hydrate)
            #add H using reduce doi.org/10.1006/jmbi.1998.2401
            #H_receptor_files = [files.replace("_clean.pdb", "_H.pdb") for files in clean_files]
        

                
        reduce=Reduce2App()
            
        results=reduce(in_files=receptors)
        print(results)
            
 
        tasks.append((self,box[_id],receptors[_id],_id,hydrate))
        
        #tasks = [(self,box[ids],files,ids,hydrate) for files,ids in zip(H_receptor_files,receptor_ids)]  #calculate task paramenters
        
        
        
        with Pool(n_proc) as pool:
            print("Start target preparation")
            results=pool.starmap(preparation, tasks)

        #sostituire questo con dataframe in shared memory
        #print(results)
        data={}
        #[{},{}] array di diz -itero array->d{},d{} se d['id'] == box['id'] 
        for column in results[0].keys():
            
            data[column] = tuple(result[column] for result in results)

        self.target_data=pd.DataFrame(data)
        self.target_data=self.target_data.set_index("id")
        print(self.target_data.head())
        
        end = time.perf_counter()
        
        print(f"Terminated in {end-start:.3f}s,avg:{(end-start)/len(results):.3f} s/receptor")


    def __dock__(self):
        pass


    def screen_single_target(self):
        pass

    def mmCIF_to_pdb(filepath:str):
        #fai quello che serve
        return Screener.pdb_fix(filepath)

    @staticmethod
    def prody_load(filepath:str):
        return parsePDB(filepath)


    """
    calculate box  based off crystallographic ligand
    """
    def cligand_box(self,filepath:str,ligand_name:str):
        base_pdb=parsePDB(filepath)
        ligand_selection=base_pdb.select(f"resname {ligand_name}") #select cristallographyc ligand for gridbox center
        #print(ligand_selection)
        cryst_lig_path=f"{self.ligands_workpath}/cryst_{ligand_name}.pdb"
        writePDB(cryst_lig_path,ligand_selection)
        lig_mol=Chem.MolFromPDBFile(cryst_lig_path,removeHs=False)
        lig_mol=Chem.AddHs(lig_mol)
        center,size=self.get_box(lig_mol)
        box={"box_src":ligand_name,"box_type":"LIG","box_center":center,"box_size":size,"box_mol":lig_mol}
        
        #center_x, center_y, center_z = calcCenter(ligand_selection) 
        return box
    
    def residue_box():
        raise NotImplementedError


    @staticmethod
    def prody_fix(filepath:str,Dir:str,Id:str,hydrate:bool=False,altloc:str="A",prody_string:str=None,*args):
        if prody_string is None:
            prody_string="chain A not hetero"
        if hydrate:
            prody_string += "and water "
        
        base_pdb=parsePDB(filepath)
        clean_pdb=f"{Dir}/{Id}/{Id}_clean.pdb"
        
        if checkNonstandardResidues(base_pdb):
            #remove non std residue (missing)
            protein=base_pdb.select(prody_string)
        else:
            protein=base_pdb.select(prody_string)
        
        #alt_prot=protein.select(f"altloc {altloc}")
        writePDB(clean_pdb,protein)
        
        
        
        return clean_pdb

    pass
    @staticmethod
    def pdb_fix(filepath:str,Dir:str,Id:str,hydrate:bool=False,*args):
        fixer = PDBFixer(filename=filepath)
        fixer.findMissingResidues() #find and repalce missing residues in the file
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues() #find and replace non-standard resisdues
        fixer.removeHeterogens(keepWater=hydrate)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        if hydrate: #prepare for hydrated docking
            fixer.addSolvent(fixer.topology.getUnitCellDimensions())
        

        clean_pdb = f"{Dir}/{Id}/{Id}_clean.pdb"  
        with open(clean_pdb, 'w') as pdb_file:
            PDBFile.writeFile(fixer.topology, fixer.positions,pdb_file)

        return clean_pdb


    @staticmethod
    def write_pdbqt(dir:str,filename:str,pdbqt):
        if pdbqt[1]:
            pdbqt=pdbqt[0]
            with open(path.join(dir, f"{filename}.pdbqt"), "w") as f:
                f.write(pdbqt)
            return pdbqt
        else:
            print(f"\nPreparation of {filename} Failed:{pdbqt[2]}")
    pass
    
    def mmCIF_to_pdb():
        raise NotImplementedError


    """
    check if molecule is in box
    """

    def _is_enveloping(mol,size:Iterable,center:Iterable=None,conf_id:int=0):
        conf = mol.GetConformer(id=conf_id)
        size = np.array(size) / 2
        coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())] #generate list with all atoms 3D coordinates
        coords = np.array(coords)
        if center is None:
            center = coords.mean(axis=0)
        min_bound = center - size 
        max_bound = center + size
        inside=np.all((coords >= min_bound) & (coords <= max_bound))
        return inside




    @staticmethod
    def get_box(mol, box_mode:str= "geometric", conf_id:int=0, padding:float=5.0): #preparation
        conf = mol.GetConformer(id=conf_id)

        coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]      #generate list with all atoms 3D coordinates

        box_calculation={"geometric":calc_box,"radius_of_gyration":radius_of_gyration}

        box=box_calculation.get(box_mode) 
        if box is None:
            raise Exception("Invalid box calculation type")
        
        
        box_dims=np.around(box(coords,padding),3)
        
        #check if output is correct
        for i, values in enumerate(box_dims, start=1):
            if not hasattr(values, '__iter__'):
                raise TypeError("Returned single box value for center or dimensions")
            if len(values) != 3:
                raise ValueError("Returned less dimension or center values then expected")
            if not all(isinstance(v, float) for v in values):
                raise TypeError("Returned bad center or dimensions value types,float is expected")
        return box_dims


    @staticmethod
    def generate_vina_config(mol,save_prefix:str,filename:str,padding:float=5.0,confId:int=0,exhaustiveness:int=24,poses:int=5)->None:
        print("Sbagliato da rivedere prima di usare")
        raise NotImplementedError
        #addictional_arguments=[]
        save_path=path.join(save_prefix,f"{filename}.txt")
        center, size = self.get_box_dimensions(mol, conf_id=confId)#Sbagliato da rivedere prima di usare
        vina_config="\n".join([f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x = {size[0]:.3f}",
        f"size_y = {size[1]:.3f}",
        f"size_z = {size[2]:.3f}",
        f"exhaustiveness = {exhaustiveness}",
        f"num_modes = {poses}"])
        with open(save_path, "w") as f:
            f.write(vina_config)
            pass
    pass



    
    def from_protein_data_bank(self,pdb_ids:Iterable[str],destination_folder:str=None,compressed:bool=False):
        paths={}
        for pdb_id in pdb_ids:
            filename = '%s.pdb' % pdb_id
            # Add .gz extension if compressed
            if compressed:
                filename = '%s.gz' % filename
            
            if destination_folder is None:
                destination_folder =  path.join(self.source_dir,"target")
            destination_file = path.join(destination_folder, filename)

            if path.exists(destination_file):
                paths[pdb_id]=destination_file
            else:
                # Download the file
                url = 'https://files.rcsb.org/download/%s' % filename
                try:
                    urlretrieve(url, destination_file)
                except Exception as e:
                    print(e)
                paths[pdb_id]=destination_file

            if compressed:
                with gzip.open(destination_file, 'rb') as f_in:
                    with open(destination_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                        pass
            pass
        return paths
        pass



    def open_sdf(filepath:str,filetype:str,base_path:str):

        if filetype == "sdf":
            file=open(filepath,"rb")
        else:
            file=gzip.open(filepath, 'rb')
        
        with Chem.ForwardSDMolSupplier(file) as supply:
            mols,ids=sdf_from_supplier(suppl=supply,basepath=base_path)    
        
        return mols,ids

    def mmCIF_to_pdb():
        pass

    def sdf_from_supplier(suppl,basepath:str):
        mols=set()
        for mol in supplier:
            if mol is None:
                continue
            
            cid=''.join(random.choice(chars) for _ in range(size)) #assign random id
            mol.SetProp("id",cid)
            ligand_path=path.join(base_path,f"{cid}")
            makedirs(ligand_path,exist_ok=True)
            outsdf_path=f"{ligand_path}/{cid}.sdf"
            mols.add(mol)
            writer = Chem.SDWriter(outsdf_path)
            writer.write(mol)
            writer.close()
        return mols


    def open_csv(filepath:str,*args):
        mols=set()
        with open(filepath, newline='\n') as file:
            reader = csv.DictReader(file)
            for row in reader:
                if row['id'].isalnum(): #da cambiare
                    ligand_name=row['id']
                else:
                    ligand_name=''.join(random.SystemRandom().choice(string.ascii_uppercase) for _ in range(6))
                    #ids.add(ligand_name) #generate random alphanumerical id
                mol=Chem.MolFromSmiles(SMILES=row['smiles'])
                if mol is None:
                    continue
                mol=Chem.AddHs(mol)
                #print(ligand_name)
                mol.SetProp("id",ligand_name)
                mols.add(mol)
        pass
        
        return mols
    
    """
    force_field:Callable force field funtion for gridbox calculation 
    hydrate: toggle hydration protocol 
    """
    def one_map_for_each(self,ff,receptor_id:str,receptor:str,hydrate:bool=False,n_conf:int=1):
        #spostare selzione algo in process 
        if ff is None:
            raise Exception(f"Missing forcefield")



        #itera attraverso i ligando con i nomi 
        # crea cartella per le mappe in ognuna di quelle dei ligandi
        #calcola mappe dei ligand types e usando la box presente nei dati dei ligandi ligands.data
        
        coupled_fn={}
        for ligand_id,ligand in self.ligand_data[["base","box_center","box_size"]].iterrows():
            save_prefix=f"{self.maps_workpath}/{ligand_id}"
            makedirs(save_prefix,exist_ok=True)
            #find best box for configuration
            #fetch ligand types
            
            coupled=ff.__map__(ligand_id,ligand["base"],receptor,receptor_id,ligand["box_center"],ligand["box_size"],save_prefix)
            coupled_fn[ligand_id] = coupled
        
        return coupled_fn


    



    def one_map_for_all():
        raise NotImplementedError
         #itera attraverso i ligando con i nomi 
        # crea cartella per le mappe unica
        #calcola mappe dei ligand types e usando la box presente nei dati dei ligandi ligands.data 
        pass

    def flexres_map():
        raise NotImplementedError
        pass

    def __set_ligands__(self,attribute,value):
        setattr(cls,attibute,value)

    def __get_ligands__(cls,attribute):
        return getattr(cls,attibute)

    def __set_target__():
        pass

    def __get_target__():
        pass

    def __set_docking_data__():
        pass

    def __get_docking_data__():
        pass

    def parse_residue_string(string:str): #From Meeko
        """
        Args:
            string (str): Residue identifier in the format 'chain:resname:resnum'.

        Returns:
            tuple: Residue identifier components: (chain, resname, resnum).

        Example:
            >>> parse_residue_string("A:ALA:42")
            (('A', 'ALA', 42), True, '')
        """

        res_id = (None, None, None)
        if string.count(":") != 2:
            raise Exception("Need exacly two ':' but found %d in '%s'" % (string.count(":"), string))
        chain, resname, resnum = string.split(":")
        if len(chain) not in (0, 1):
            raise Exception("chain must be 0 or 1 char but it is '%s' (%d chars) in '%s'" % (chain, len(chain), string))
        if len(resname) > 3:
            raise Exception("resname must be max 3 characters long, but is '%s' in '%s'" % (resname, string)) 
        try:
            resnum = int(resnum)
        except:
            raise Exception("resnum could not be converted to integer, it was '%s' in '%s'" % (resnum, string)) 
        res_id = (chain, resname, resnum)
        return res_id
    
    def no_map(self,receptor_id:str):
        
        pass

    """
    screeen_single_target:
        doking_agent:Callable docking function and saves results
        map_calculator:callable map calculator
    """


    def coords_align(mol,align_pivot:Iterable,conf_id:int=-1):
        
        mol_conf=mol.GetConformer(conf_id)
        coords = np.array([list(mol_conf.GetAtomPosition(i)) for i in range(mol_conf.GetNumAtoms())])
        center,_=Screener.get_box(mol,conf_id=conf_id)
        T=align_pivot-center #get translation vector
        translated_coords=coords+T #apply translation vector
        for i in range(mol.GetNumAtoms()):
            x,y,z = translated_coords[i]
            mol_conf.SetAtomPosition(i,Point3D(x,y,z))
        return mol  


    def analyze(self,single_target:bool=True):
        if not single_target: #gestione multi target
            pass
        












pass



smiles=[
    "HU-210", "CCCCCCC(C)(C)C1=CC2=C([C@@H]3CC(CO)=CC[C@H]3C(C)(C)O2)C(C)=C1",
    "(−)-11-nor-9-Carboxy-Δ9-THC","[H][C@@]12CC(=CC[C@@]1([H])C(C)(C)OC1=C2C(O)=CC(CCCCC)=C1)C(O)=O",
    "(+/-)-11-nor-9-Carboxy-Δ9-THC","[2H]C([2H])([2H])CCCCC1=CC(=C2[C@@H]3C=C(CC[C@H]3C(OC2=C1)(C([2H])([2H])[2H])C([2H])([2H])[2H])C(=O)O)O",
    "DHT","C[C@]12CCC(=O)C[C@@H]1CC[C@@H]3[C@@H]2CC[C@]4([C@]3(CC[C@@H]4O)C)C",
    "Testosterone","C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C",
    "CBD","CCCCCC1=CC(O)=C([C@@H]2C=C(C)CC[C@H]2C(C)=C)C(O)=C1",
    "CBN","CCCCCC1=CC2=C(C(O)=C1)C1=CC(C)=CC=C1C(C)(C)O2",
    "R1881","C[C@@]1(CC[C@@H]2[C@@]1(C=CC3=C4CCC(=O)C=C4CC[C@@H]23)C)O",
    "THC","CCCCCC1=CC(=C2[C@@H]3C=C(CC[C@H]3C(OC2=C1)(C)C)C)O",
    "5F-AKB48","C1C2CC3CC1CC(C2)(C3)NC(=O)C4=NN(C5=CC=CC=C54)CCCCCF",
    "MDMB-CHMINACA","CC(C)(C)[C@@H](C(=O)OC)NC(=O)C1=NN(C2=CC=CC=C21)CC3CCCCC3"
]



#if __name__ == "__main__":




#if __name__ == "__main__":
#agent per il calcolo delle mappe
#agent per il docking
#screen calcola le mappe nel protocollo corretto e docka







"""
The majority of AR mutations result in single amino acid
substitutions, which are mostly found in the AR androgen-binding domain. The mutation
T877A, which has been found in roughly 30% of metastatic CRPC patients, is the most
common [134]. Other mutations have resulted in enhanced AR binding to coregulators,
resulting in higher AR transcriptional activity vis-à-vis H874Y and W435L mutations.
These mutations have been implicated in the development of AR resistance arising from AR-targeted therapy [124].


"""