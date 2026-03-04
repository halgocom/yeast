from abc import ABC,abstractmethod

class ScreenTool(ABC):
    require_maps:bool=True

    @abstractmethod
    def __dock__(self,*args,**kwargs):
        ...
    
    @abstractmethod
    def __map__(self,*args,**kwargs):
        ...