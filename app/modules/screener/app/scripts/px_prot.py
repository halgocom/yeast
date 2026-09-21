from filters import LipinskiRDKit
from forcefield import RdkitEtkgd
from protenix_pipeline import ProtenixPipeline
from ProtenixScreen import calc_maps,pxdock

if __name__ == "__main__":
    clean_settings = {"5s8i_ligand":None,"5s8i_chain":"A","5s8i_altloc":"A"}#,"pH":6.1,"forcefield":"amber"}

    dock=ProtenixPipeline(experiment_id="prost_canc",replicates=5,run_name="protenix")
    
    dock.load_targets(from_PDB=["5s8i"])
    #dock.target_cleanup(**clean_settings,repair=False)
    dock.process_targets(box={"5s8i":((27.333,1.472,3.639),(24,32,16))})

    dock.load_ligands(filename="cannabinoids.csv",sanitize=False)
    
    dock.filter_ligands(filter_f=LipinskiRDKit,allow_partial=True)
    
    dock.optimize_ligand(opt_f=RdkitEtkgd,n_conf=5,pH_max=6.5,pH_min=6.1)
    
    dock.process_ligands(output="path")
    
    dock.init_docking()
    #dock.set_maps(map_f=calc_maps,map_mode="o4e",target_id="2am9",target_box=True)
    #print(dock.map_data)
    dock.screen_single_target(dock_f=pxdock,target_id="5s8i")
    