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
from os import makedirs,getcwd

def prepare_ligand(mol,out_dir,ligand_id:str,conf:int=0,hydrate:bool=False,save_sdf:bool=True,preparation_config:dict=None,*args,**kwargs):
    #mol = Chem.Mol(mol)

    
    # Molecule preparation from Meeko
    if preparation_config is not None:
        preparator = MoleculePreparation.from_config(preparation_config) #load config
    else:   
        preparator = MoleculePreparation(hydrate=hydrate) #hydrate protocol

    #makedirs(ligand_dir,exist_ok=True)

    #Chem.rdMolAlign.AlignMolConformers(mol)
    
    

    sizes=[]
    confs=[]
    sdf=[]
    for conf_id in range(mol.GetNumConformers()):
        
        center,conf_size=get_box(mol=mol,conf_id=conf_id)
        
        
         
        conf=preparator.prepare(mol,conformer_id=conf_id)
        
        conf_pdbqt = PDBQTWriterLegacy.write_string(conf[0])
        write_pdbqt(dir=out_dir,filename=f"{ligand_id}_{conf_id}",pdbqt=conf_pdbqt) #save conf pdbqt
        
        
        
        
        confs.append(f"{out_dir}/{ligand_id}_{conf_id}.pdbqt")

            

    


    center = [float(ele) for ele in center]
    size = [float(ele) for ele in conf_size]
    return {"ligand_id":f"{ligand_id}","conf_pdbqt":tuple(confs),"center":center,"size":size}




ligs=["inib.sdf","r1881.sdf"]
ligands = {}
for lig in ligs:
    lig_path = path.join(getcwd(),"lig_prep",lig)
    mol = Chem.SDMolSupplier(lig_path)[0]
    mol=Chem.AddHs(mol,addCoords=True)
    out_dir = path.join(getcwd(),"lig_prep")
    _id = lig.split(".")[0]

    preped=prepare_ligand(mol,out_dir,_id)
    json_file=path.join(getcwd(),"lig_prep",f"{_id}.json")
    print(json_file)
    with open(json_file,"w") as f:
        json.dump(preped,f)