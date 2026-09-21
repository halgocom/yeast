from abc import ABC,abstractmethod
from dataclasses import dataclass
from typing import Iterable

class EnvApp(ABC):


    """
        EnvApp.__init__:initializes class with information about used env
        env_name:str name of used enviroment
    """

    @abstractmethod
    def __init__(self,env_name:str):
        ...
    """
        EnvApp.__call__:calls and execute the scripts in an external conda env
            script_path:str absolute path to script to execute
            params: list of script paramenter for subprocess call
            concurrent:bool toggle for concurrent execution
    """
    @abstractmethod
    def __call__(self,*args,**kwargs):
        ...

class CondaExec(EnvApp):

    @abstractmethod
    def __init__(self,env_name:str):
        ...

    @abstractmethod
    def __call__(self):
        ...

    #method to grab env python interpreter path
    def __get_py__(env_name:str,env_type:str):
        pass

class AppInterface(ABC):
    
    app_path = "..."
    