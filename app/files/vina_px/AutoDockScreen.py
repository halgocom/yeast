from vina import Vina
from typing import Iterable

import os
from numpy import round
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate


def calc_map(target:str,
    target_id:str,
    box_center:Iterable,
    box_size:Iterable, 
    save_prefix:str,
    ligand_id:str=None,
    sf_name='vina',
    cpu:int=0,
    seed=0,
    no_refine:bool=False,
    verbosity:int=1,**kwargs):

    os.makedirs(save_prefix, exist_ok=True)  

    exe = Vina(sf_name=sf_name,
               cpu=cpu,
               seed=seed,
               no_refine=no_refine,
               verbosity=verbosity)

    

    map_log = {
        "target_id": target_id,
        "box_center": box_center,
        "box_size": box_size
    }

    
    exe.set_receptor(target)

    exe.compute_vina_maps(list(map(float, box_center)),
    list(map(float, box_size)),
    force_even_voxels=True)

    if ligand_id:
        fn = f"{target_id}_{ligand_id}"
        map_log["ligand_id"]=ligand_id
    else:
        fn = target_id
        

    coupled_fn = f"{save_prefix}/{fn}"

    exe.write_maps(coupled_fn, overwrite=True,fld_filename=fn)
    map_log["coupled_fn"] = coupled_fn



    return map_log

def dock(exp_type:str="basic",**kwargs):
    exp_type=exp_type.lower()
    docking_modes={"basic":__basic_docking__,"hydrate":__hydro_docking__,"flexres":__flex__}
    dock=docking_modes.get(exp_type)
    if dock is None:
        supported=docking_modes.keys().join(",")
        raise RuntimeError(f"Unsupported docking mode:{exp_type},\nSupported modes {supported}")
    return dock(**kwargs)


def basic_docking(
    ligand_id:str,
    ligand_paths:str,
    target_id:str,
    coupled_fn:str,
    save_path:str,
    replicate_id:int,
    conf_id:int=0,
    poses:int=20,
    write_poses:int=1,
    sf_name='vina',
    cpu:int=0,
    seed=124657,
    no_refine:bool=False,
    verbosity:int=0,
    exhaustiveness:int=32,
    ):
    
    
    save_path=f"{save_path}/{ligand_id}"
    os.makedirs(save_path, exist_ok=True) 
    exe = Vina(sf_name=sf_name,
        cpu=cpu,
        seed=seed,
        no_refine=no_refine,
        verbosity=verbosity)
    
    exe.load_maps(coupled_fn)
    
    dock_log ={"replicate":replicate_id,"ligand_id":ligand_id,"target_id":target_id,"seed":seed}



    convertion_fails = 0
    
        
    ext_filename=f"{target_id}_{ligand_id}_{conf_id}_replicate{replicate_id}"

    ext_path=f"{save_path}/{ext_filename}"
    

    exe.set_ligand_from_file(ligand_paths)
    energy = exe.score()


    # Minimized locally the current pose
    energy_minimized = exe.optimize()


    exe.write_pose(f'{ext_path}_min.pdbqt', overwrite=True)

    # Dock the ligand
    exe.dock(exhaustiveness=exhaustiveness, n_poses=poses)
    exe.write_poses(f'{ext_path}.pdbqt', n_poses=write_poses, overwrite=True)
    energies = exe.energies(n_poses=write_poses)

    pdbqt_mol = PDBQTMolecule.from_file(f'{ext_path}.pdbqt',skip_typing=True)


    sdf_string,failure = RDKitMolCreate.write_sd_string(pdbqt_mol)


    if sdf_string == "" or len(failure) !=0:
        convertion_fails+=1
    else:
        with open(f'{ext_path}.sdf',"w") as f:
            f.write(sdf_string)

        docked_pdbqt=f'{ext_path}.pdbqt'
        bind_energies=energies[0][0]
        docked_sdf=f'{ext_path}.sdf'

    dock_log.update({"conf_id":conf_id,f"docked_pdbqt":docked_pdbqt,f"docked_sdf":docked_sdf,f"energies":bind_energies,"fails":convertion_fails})

    return dock_log
    

def __hydro_docking__(*args,**kwargs):
    raise NotImplementedError
    

def __flex__(*args,**kwargs):
    raise NotImplementedError






"""
def calc_map(receptor:str,receptor_id:str,box_center:Iterable,box_size:Iterable,save_prefix:str,ligand_id:str=None,redock_box_delta:Iterable=None,redocking:bool=False,sf_name='vina',cpu:int=0,seed=0,no_refine:bool=False,verbosity:int=0,**kwargs):
        exe=Vina(sf_name=sf_name,cpu=cpu,seed=seed,no_refine=no_refine,verbosity=verbosity)
        
        map_log ={"target_id":receptor_id,"box_center":box_center,"box_size":box_size}
        
        redock_box_size =box_size+redock_box_delta
        
        exe.set_receptor(receptor)
        exe.compute_vina_maps(box_center,box_size,force_even_voxels=True)
        if not ligand_id:
            fn =f"{receptor_id}"
            coupled_fn=f"{save_prefix}/{receptor_id}"
            exe.write_maps(coupled_fn,overwrite=True,fld_filename=fn)
            map_log =map_log | {"coupled_fn":coupled_fn}
        else:
            fn = f"{receptor_id}_{ligand_id}"
            coupled_fn=f"{save_prefix}/{ligand_id}"
            exe.write_maps(coupled_fn,overwrite=True,fld_filename=fn)
            map_log =map_log | {"coupled_fn":coupled_fn}
        
        if redocking:
            redock_fn=f"{fn}_redock"
            redock_coupled_fn=f"{coupled_fn}_redock"
            exe.compute_vina_maps(box_center,redock_box_size,force_even_voxels=True)
            exe.write_maps(redock_coupled_fn,overwrite=True,fld_filename=redock_fn)
            map_log = map_log | {"redock_coupled_fn":redock_coupled_fn,"redock_box_size":redock_box_size}
        return map_log







class VinaScreen(ScreenTool):
    def __init__(self,sf_name='vina',cpu:int=0,seed=0,no_refine:bool=False,verbosity:int=0,exp_type:str="basic",*args,**kwargs):
        
        
        self.init_kwargs = {"sf_name":sf_name,"cpu":cpu,"seed":seed,"no_refine":no_refine,"verbosity":verbosity}


    maps:calculate map of single target
    ligand_id:ligand_id str
    ligand:ligand pdbqt str
    receptor:target pdbqt filepath str
    receptor_id:target id str

    *args
    

    
    


    def __get_exe__(self):
        return self.EXE

    def reset(self,**kwargs):
        self.init_kwargs.update(kwargs)
        self.EXE=Vina(**self.init_kwargs)

"""