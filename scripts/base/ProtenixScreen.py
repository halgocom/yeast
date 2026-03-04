from docking import ScreenTool
from interface import CondaExec
from typing import Iterable
from dotenv import load_dotenv
from os import getenv
import subprocess
import os



def calc_maps(
    target:str,
    target_id:str,
    box_center:Iterable,
    box_size:Iterable, 
    save_prefix:str,
    ligand_id:str=None,
    threads:int=0):
    
    os.makedirs(save_prefix, exist_ok=True)  


    

    map_log = {
        "target_id": target_id,
        "box_center": list(map(float, box_center)),
        "box_size": list(map(float, box_size))
    }


    if ligand_id:
        fn = f"{target_id}_{ligand_id}"
        map_log["ligand_id"]=ligand_id
    else:
        fn = target_id
        

    coupled_fn = f"{save_prefix}/{fn}"

    map_log["coupled_fn"] = coupled_fn



    return map_log


def pxdock(receptor_id:str,receptor:str,ligand:str,ligand_id:str,replicate_id:int,save_path:str,box_center:Iterable,box_size:Iterable,conf_id:int,seed:int,threads:int,):

    
    box_center = [str(coord) for coord in box_center]

    box_size = [str(coord) for coord in box_size]

    seed = str(seed)
    box_center = ",".join(box_center)
    box_size=",".join(box_size)

    replicate_id = str(replicate_id)
    conf_id = str(conf_id)
    threads = str(threads)
    proc=subprocess.run(
        ["micromamba", "run", "-n", "pxdock", "python", "/home/screener/modules/pxdock/screen.py",receptor_id, receptor, ligand,ligand_id,replicate_id,save_path,box_center,box_size,conf_id,seed,threads],
        
        capture_output = True, 
        text = True 
    )
    print(proc.stdout)
    print(proc.stderr)

    return "1"













"""
class ProtenixScreen(ScreenTool,CondaExec):


    def __init__(self,env_name:str="pxdock",env_type:str="micromamba"):
        self.env_type=env_type
        load_dotenv(dotenv_path="/home/screener/env/.env")
        module_path=getenv("MODULE_PATH")
        self.env_name=env_name 
        self.script_path = f"{module_path}/{env_name}/screen.py"
        self.box_center={}
        self.box_size={}
        self.maps_prefix=None
        self.ligands={}
        pass


    def __dock__(self,ligand_pdbqt:str,ligand_id:str,receptor_pdbqt:str,receptor_id:str,replicate_id:int,save_path:str,coupled_fn:str,conf_id:int=-1):
        conf_id = str(conf_id)
        replicate_id = str(replicate_id)
        return ProtenixScreen.__call__(self,receptor_id=receptor_id,receptor=receptor_pdbqt,ligand=ligand_pdbqt,ligand_id=ligand_id,coupled_fn=coupled_fn,conf_id=conf_id,replicate_id=replicate_id,save_path=save_path)
    



    #method to grab env python interpreter path
    def __get_py__(env_name:str,env_type:str):
        pass

    def reset_seed(self):
        ...
"""
if __name__=="__main__":
    main()