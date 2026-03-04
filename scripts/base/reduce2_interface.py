from interface import CondaExec
import subprocess
from dotenv import load_dotenv
from typing import Iterable
from os import getenv
from multiprocessing import pool
from pathlib import Path


def addH(env_name:str,env_type:str,script_path:str,receptor_id:str,in_file:str):
    
    dirname = Path(in_file).parent
    out_file = f"{dirname}/{receptor_id}_H.pdb"
    print(f"Adding H to receptor {receptor_id}")
    try:
        subprocess.run(
            [env_type, "run", "-n", env_name, "python", script_path, in_file, out_file],
            check=True
        )
    except Exception as e:
        print(f"Adding H failed for {receptor_id} \n")
        print(e)
        return None
    return {f"{receptor_id}":out_file}





class Reduce2App(CondaExec):

    def __init__(self,env_name:str="CCTBX",env_type:str="micromamba"):
        #get necessary info
        self.env_type=env_type
        load_dotenv(dotenv_path="/home/screener/env/.env")
        module_path=getenv("MODULE_PATH")
        self.env_name=env_name 
        self.script_path = f"{module_path}/{env_name}/run_reduce.py"
        
    
    def __call__(self,in_files:dict,*args):
        
        rec_h={}
 
        for _id in in_files.keys():
            rec_h[_id]=addH(self.env_name,self.env_type,self.script_path,_id,in_files[_id])

        return rec_h


    
#if __name__ == "__main__":
    #red=Reduce2App()
    #infles={"2am9":"/home/screener/files/prost_canc/targets/2am9/2am9_clean.pdb"}
    #result=red(in_files=infles,concurrent=True)
    #print(result)

pass




pass