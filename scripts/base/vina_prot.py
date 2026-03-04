from mk_pipeline import MeekoPipeline
from filters import LipinskiRDKit
from forcefield import RdkitEtkgd
from AutoDockScreen import calc_map,basic_docking
import os




if __name__ == "__main__":

        


    clean_settings = {"2am9_ligand":"TES","2am9_chain":"A","2am9_altloc":"A","pH":6.1,"forcefield":"amber"}

    dock=MeekoPipeline(experiment_id="prost_canc",replicates=1,run_name="time_comp")

    
    
    dock.load_targets(from_PDB=["2am9"])
    dock.target_cleanup(**clean_settings)
    dock.process_targets(box={"2am9":((27.333,1.472,3.639),(24,32,16))})


    
    dock.load_ligands(filename="cannabinoids.csv",sanitize=False)
    
    dock.filter_ligands(filter_f=LipinskiRDKit,allow_partial=True)
    
    dock.optimize_ligand(opt_f=RdkitEtkgd,n_conf=5,pH_max=6.5,pH_min=6.1)
    
    dock.process_ligands(output="path")
    
    
    dock.init_docking()

    
    dock.calculate_maps(map_f=calc_map,map_mode="o4a",target_id="2am9",target_box=True)

    dock.screen_single_target(dock_f=basic_docking,target_id="2am9",verbosity=0)

    #dock.check_poses(target_id="2am9",redock="Testosterone")
    #dock.analyze_interactions(target_id="2am9")
    #dock.calc_metrics(("ki","le"),target_id="2am9")

    #dock.create_report
