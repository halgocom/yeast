from mk_pipeline import MeekoPipeline
from filters import LipinskiRDKit
from forcefield import RdkitEtkgd
from AutoDockScreen import calc_map,basic_docking





if __name__ == "__main__":

    
    clean_settings = {"1z95":{"ligand":"198","chain":"A","altloc":"A","pH":6.1}}

    dock=MeekoPipeline(experiment_id="prost_canc",replicates=5,run_name="vina_inib")

    
    
    dock.load_targets("1z95")
    dock.target_cleanup(**clean_settings)
    dock.process_targets({"1z95": ((28.1117, 1.4037, 6.5190),(26.0, 26.0, 26.0))})

    
    dock.load_ligands(filename="cannabinoids_inib_opt.sdf.gz",sdf_confs=True)
    
    dock.filter_ligands(filter_f=LipinskiRDKit,allow_partial=True)
    
    dock.optimize_ligand(protonate=True,pH_max=6.5,pH_min=6.1)

    dock.process_ligands(output="path")
    
    dock.init_docking(safety_factor=1)

    dock.calculate_maps(map_f=calc_map,map_mode="o4a",target_id="1z95",target_box=True)

    dock.screen_single_target(dock_f=basic_docking,target_id="1z95",verbosity=0,sf_name="vinardo",poses=10)

    dock.check_poses(target_id="1z95",redock="RBICALUTAMIDE",box_threshold=0.9)
    dock.analyze_interactions(target_id="1z95")
    dock.calc_metrics(("ki","le"),target_id="1z95")

    dock.create_report