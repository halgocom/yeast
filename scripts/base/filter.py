
from abc import ABC,abstractmethod
from dataclasses import dataclass


class Filter(ABC):
    def __init__(self,*args,**kwargs):
        ...
    
    @abstractmethod
    def __filter__(self,*args,**kwargs):
        ...
