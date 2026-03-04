from pxdock import ProtenixDock
import sys
import json
from pathlib import Path
import os
"""
def save(json:dict,save_path:str):
    with open(save_path,"w") as f:
        json.dump(json,f)
"""

def main():
    
    print(sys.argv)
    receptor_id = sys.argv[1]
    receptor = sys.argv[2]
    ligand = sys.argv[3]
    ligand_id = sys.argv[4]
    replicate_id = sys.argv[5]
    save_path = sys.argv[6]
    in_box_center =sys.argv[7]
    in_box_size = sys.argv[8]
    conf_id = sys.argv[9]
    seed = sys.argv[10]
    threads = sys.argv[11]
    
    
 #receptor_id:str,receptor:str,ligand:str,ligand_id:str,replicate_id:int,save_path:str,box_center:Iterable,box_size:Iterable,conf_id:int,threads:int
    seed = int(seed)
    save_lig_path = os.path.join(save_path,ligand_id)
    os.makedirs(save_lig_path,exist_ok=True)
    box_center = [float(value) for value in in_box_center.split(",")]
    box_size = [float(value) for value in in_box_size.split(",")]#convert back box center and box size
    
    #box_center = [0., 0., 0.] # box center for receptor
    #box_size = [10., 10., 10.] # box size for receptor
    dock_instance = ProtenixDock(receptor,nthreads=threads)
    dock_instance.set_box(box_center, box_size)

    # Optional: you can generate cache maps for receptor, and then you can load it for next docking.
    # In our tests, setting this parameter to 0.175 can achieve a balance between effect and performance.
    # out_dir = dock_instance.generate_cache_maps(spacing=0.175)
    # and in next run:
    # dock_instance.load_cache_maps(out_dir)

    # the docking_res_files is in json format.
    docking_res_files = dock_instance.run_docking(ligand,out_dir=save_lig_path)
    #save(docking_res_files,f"{save_path}/{receptor_id}_{ligand_id}_replicate{replicate_id}.json")

    print(docking_res_files)
    
    
    pass

if __name__ == "__main__":
    main()