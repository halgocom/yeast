from AutoDockScreen import *

from os import makedirs,getcwd,path
import json
jobs = {

    "inib":{
        "ligand":"/home/screener/files/vina_px/docking/inib_0.pdbqt",
        "target":"/home/screener/files/vina_px/docking/1z95_inib.pdbqt",
        "box":[[-5.238, 2.604, 3.253],[19.075, 22.91, 16.975]]
    },
   "r1881":{
        "ligand":"/home/screener/files/vina_px/docking/r1881_0.pdbqt",
        "target":"/home/screener/files/vina_px/docking/1z95_r1881.pdbqt",
        "box": [[-4.889, -6.0, 0.012],[15.209, 20.878, 17.309]]
    }
}

base_dir = path.join(getcwd(),"docking")
for i in range(10):

    for job in jobs.keys():

        job_path = path.join(base_dir,job)
        makedirs(job_path,exist_ok=True)
        map_path = path.join(job_path,"maps")
        makedirs(map_path,exist_ok=True)
        replicate_path =path.join(job_path,f"{i}")
        
        exe = Vina(sf_name="vina",
            cpu=0,
            seed=0)

        exe.set_receptor(jobs[job]["target"])
        exe.compute_vina_maps(jobs[job]["box"][0],jobs[job]["box"][1],force_even_voxels=True)
        coupled_fn = f"{map_path}/{job}"
        exe.write_maps(coupled_fn, overwrite=True)

        exe.set_ligand_from_file(jobs[job]["ligand"])
        energy = exe.score()

        exe.dock(exhaustiveness=32,n_poses=20)

        dock_path = path.join(job_path,"dock",job)
        makedirs(dock_path,exist_ok=True)
        ext_filename = path.join(dock_path,f"{job}_{i}.pdbqt")
        exe.write_poses(ext_filename, n_poses=1, overwrite=True)
        energies = exe.energies(n_poses=1)

        pdbqt_mol = PDBQTMolecule.from_file(ext_filename,skip_typing=True)


        sdf_string,failure = RDKitMolCreate.write_sd_string(pdbqt_mol)

        docked_sdf=ext_filename.replace("pdbqt","sdf")

        with open(docked_sdf,"w") as f:
            f.write(sdf_string)

        
        bind_energies=energies[0][0]

        jobs[job][f"replicate {i}"] = {"sdf":docked_sdf,"pdbqt":ext_filename,"energy":bind_energies}

print(jobs)


json_file=path.join(job_path,"result.json")

with open(json_file,"w") as f:
    json.dump(jobs,f)
        




        
    

