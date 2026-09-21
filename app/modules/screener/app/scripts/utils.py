from os import path
import subprocess
import sys


#defines 
class ExpData():
    def __init__():
        pass


class Analize():
    pass


class Script():
    def __init__(self,in_filedir:str,out_filedir:str,script_dir:str,script_name:str):
        self.in_file_dir=in_file_dir
        self.out_filedir=out_filedir
        self.script_name=f"{script_name}.py"
        self.python_path=str(sys.executable) #python executable path, for _script

        self.script_path=f"{path.join(script_dir,self.script_name)}"
    
    def __call__(self,**kwargs):
        print(kwargs)
        print(self.python_path,self.script_path)
        cmd=[f"{self.python_path}",self.script_path]
        #if parameters is not None:
            #cmd.extend(parameters)
        result=subprocess.run(cmd,capture_output=True)
        print(result)
        pass

scr=Script("/home/screener/app/scripts/usr","mk_export")

scr(i=2,k=2)