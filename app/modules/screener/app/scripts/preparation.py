from abc import ABC,abstractmethod
from dataclasses import dataclass


@dataclass
class LigandTemplate(ABC):
    @abstractmethod
    def __init__(self,*args,**kwargs):
        ...

    @abstractmethod
    def __next__(self,*args,**kwargs):
        ...

    @abstractmethod
    def __iter__(self,*args,**kwargs):
        ...


@dataclass
class TargetTemplate(ABC):
    @abstractmethod
    def __init__(self,*args,**kwargs):
        ...

    @abstractmethod
    def __next__(self,*args,**kwargs):
        ...
        
    @abstractmethod
    def __iter__(self,*args,**kwargs):
        ...




class PreparationTemplate(ABC):

    @abstractmethod
    def process_targets(self,*args,**kwargs):
        ...
    
    @abstractmethod
    def process_targets(self,*args,**kwargs):
        ...

    @abstractmethod
    def __ligand__(self,*args,**kwargs):
        ...

    @abstractmethod
    def __target__(self,*args,**kwargs):
        ...





    



