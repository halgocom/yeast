from mk_pipeline import MeekoPipeline
from filters import LipinskiRDKit
from forcefield import RdkitEtkgd
from AutoDockScreen import calc_map,basic_docking
import os




if __name__ == "__main__":

        


    clean_settings = {"1z95_ligand":"198","1z95_chain":"A","1z95_altloc":"A","pH":6.1}

    dock=MeekoPipeline(experiment_id="prost_canc",replicates=5,run_name="vina_inib")

    
    
    dock.load_targets(from_PDB=["1z95"])
    dock.target_cleanup(**clean_settings)
    dock.process_targets(box={"1z95":((29.22,3.556,8.139),(40,34,28))})


    
    dock.load_ligands(filename="cannabinoids_inib.csv",sanitize=False)
    
    dock.filter_ligands(filter_f=LipinskiRDKit,allow_partial=True)
    
    dock.optimize_ligand(opt_f=RdkitEtkgd,n_conf=5,pH_max=6.5,pH_min=6.1)
    
    dock.process_ligands(output="path")
    
    
    dock.init_docking()

    
    dock.calculate_maps(map_f=calc_map,map_mode="o4a",target_id="1z95",target_box=True)

    dock.screen_single_target(dock_f=basic_docking,target_id="1z95",verbosity=0)

    dock.check_poses(target_id="1z95",redock="RBICALUTAMIDE")
    dock.analyze_interactions(target_id="1z95")
    dock.calc_metrics(("ki","le"),target_id="1z95")

    dock.create_report