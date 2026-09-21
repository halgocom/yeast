from posebusters import PoseBusters
from rdkit.Chem import MolFromPDBFile
from rdkit.Chem import AllChem,SDMolSupplier,SanitizeMol
import pandas as pd
from os import getcwd,path
import numpy as np
import json





jobs ={"inib":{
    "sdf":["/home/screener/files/vina_px/docking/testo/dock/testo/testo_0.sdf","/home/screener/files/vina_px/docking/testo/dock/testo/testo_1.sdf","/home/screener/files/vina_px/docking/testo/dock/testo/testo_2.sdf","/home/screener/files/vina_px/docking/testo/dock/testo/testo_3.sdf","/home/screener/files/vina_px/docking/testo/dock/testo/testo_4.sdf","/home/screener/files/vina_px/docking/testo/dock/testo/testo_5.sdf","/home/screener/files/vina_px/docking/testo/dock/testo/testo_6.sdf","/home/screener/files/vina_px/docking/testo/dock/testo/testo_7.sdf","/home/screener/files/vina_px/docking/testo/dock/testo/testo_8.sdf","/home/screener/files/vina_px/docking/testo/dock/testo/testo_9.sdf"],
    "energy":[-11.091,-11.117,-11.092,-11.096,-11.117,-11.114,-11.115,-11.116,-11.091,-11.088],
    "prot":"/home/screener/files/vina_px/rec_prep/2am9_testo.pdb"
    },
    "r1881":{
        "sdf":["/home/screener/files/vina_px/docking/cbn/dock/cbn/cbn_0.sdf","/home/screener/files/vina_px/docking/cbn/dock/cbn/cbn_1.sdf","/home/screener/files/vina_px/docking/cbn/dock/cbn/cbn_2.sdf","/home/screener/files/vina_px/docking/cbn/dock/cbn/cbn_3.sdf","/home/screener/files/vina_px/docking/cbn/dock/cbn/cbn_4.sdf","/home/screener/files/vina_px/docking/cbn/dock/cbn/cbn_5.sdf","/home/screener/files/vina_px/docking/cbn/dock/cbn/cbn_6.sdf","/home/screener/files/vina_px/docking/cbn/dock/cbn/cbn_7.sdf","/home/screener/files/vina_px/docking/cbn/dock/cbn/cbn_8.sdf","/home/screener/files/vina_px/docking/cbn/dock/cbn/cbn_9.sdf"],
        "energy":[-7.681,-7.611,-7.72,-7.682,-7.657,-7.715,-7.682,-7.643,-7.659,-7.677],
        "prot":"/home/screener/files/vina_px/rec_prep/2am9_cbn.pdb",

    }
}

check_columns=[
'mol_pred_loaded', 'mol_cond_loaded', 'sanitization',
'inchi_convertible', 'all_atoms_connected', 'no_radicals',
'bond_lengths', 'bond_angles', 'internal_steric_clash',
'aromatic_ring_flatness', 'non-aromatic_ring_non-flatness',
'double_bond_flatness', 'internal_energy',
'protein-ligand_maximum_distance', 'minimum_distance_to_protein','volume_overlap_with_protein']

result={}

for job in jobs.keys():
    clash_checker = PoseBusters(config="dock")
    
    target_mol = MolFromPDBFile(jobs[job]["prot"],sanitize=False,removeHs=False)
    docked_clash = clash_checker.bust(mol_pred=jobs[job]["sdf"],mol_cond=target_mol)
    docked_clash["energy"] = pd.Series(jobs[job]["energy"], index=docked_clash.index)
    

    dock_mask = np.all(docked_clash[check_columns].to_numpy(dtype=bool), axis=1)
    
    best_confs = docked_clash.loc[dock_mask]#, ["file","energy"]]

    file = best_confs["energy"].idxmin()
    energy = best_confs.loc[file, "energy"]

    result[job] ={"file":file,"energy":energy}

print(result)
docked_file = result["testo"]["file"][0]
original = "/home/screener/files/vina_px/check/2am9_ligand.sdf"

cryst_mol = SDMolSupplier(original)[0]
docked_mol = SDMolSupplier(docked_file)[0]

SanitizeMol(cryst_mol)
SanitizeMol(docked_mol)

o3a = AllChem.GetO3A(docked_mol, cryst_mol)
rmsd = o3a.Align()


result["testo"]["redock RMSD"]= rmsd
print(rmsd)
file_path = path.join(getcwd(),"check")
json_file=path.join(file_path,"result.json")

with open(json_file,"w") as f:
    json.dump(jobs,f)
        