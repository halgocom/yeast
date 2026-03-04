from typing import Callable,Iterable
def one_map_for_each(map_func:Callable,ligand_id:str,receptor_id:str,receptor:str,out_path:str,box:Iterable,redocking:bool,dim_delta:Iterable=None,center_delta:Iterable=None,*maps_args,**maps_kwargs):
     
    coupled_fn={}
    for ligand_id,ligand in self.ligand_data[["base","box_center","box_size"]].iterrows():
        save_prefix=f"{self.maps_workpath}/{ligand_id}"
        makedirs(save_prefix,exist_ok=True)
        coupled=map_func(ligand_id,receptor,receptor_id,box["box_center"],box["box_size"],save_prefix,*maps_args,**map_kwargs)
        
    
    return {"ligand_id":ligand_id,"receptor_id":receptor_id,"coupled_fn":coupled}


def one_map_for_all(map_func:Callable,ligand_id:str,receptor_id:str,receptor:str,out_path:str,box:Iterable=None,*maps_args,**maps_kwargs):
    if not box:
        raise RuntimeError("Missing box")

    save_prefix=f"{out_path}/{receptor_id}"
    makedirs(save_prefix,exist_ok=True)
    coupled=map_func(ligand_id,receptor,receptor_id,box["box_center"],box["box_size"],save_prefix,*maps_args,**map_kwargs)
    
    return {"ligand_id":ligand_id,"receptor_id":receptor_id,"coupled_fn":coupled}