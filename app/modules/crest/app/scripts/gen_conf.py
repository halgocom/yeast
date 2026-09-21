import sys
import os
from utils import read_csv,optimize_wrapper,read_xyz_n,crest_wrapper
from pathlib import Path
import shutil
import multiprocessing as mp
import gzip
from rdkit import Chem
from tqdm import tqdm
from functools import partial


FILE_PATH = Path(os.environ.get("SHARE_FILE_PATH"))
ORIG_FILE = FILE_PATH / f"sources/ligands/{sys.argv[1]}.csv"
#EXP = sys.argv[2]
#RUN = sys.argv[3] 
ALGO = "gfn2" 
WINDOW = 5.0
N_CONF= 10

OUT_XYZ = FILE_PATH /  f"sources/ligands/{sys.argv[1]}_opt.xyz"
OUT_SDF = FILE_PATH /  f"sources/ligands/{sys.argv[1]}_opt.sdf.gz"

n_proc = mp.cpu_count() -2

results=[]



#optizime single molecule with xtb in parallel
with open(OUT_XYZ, "w") as f, mp.Pool(processes=n_proc) as pool:

    for result in tqdm(pool.imap_unordered(optimize_wrapper,read_csv(ORIG_FILE)),desc="Minimizing ligands"):
        if result:
            results.append(result)
            f.write(result["xyz"] + "\n")
            f.flush()
pass

crest_dir = FILE_PATH / f"sources/ligands/{sys.argv[1]}"
os.makedirs(crest_dir, exist_ok=True)

crest = partial(
    crest_wrapper,
    crest_dir=crest_dir,
    algo=ALGO,
    ewin=WINDOW,
    n_conf=N_CONF,
    threads=2
)

n_proc = 2
with gzip.open(OUT_SDF, "wt") as gz:
    writer = Chem.SDWriter(gz)

    with mp.Pool(processes=n_proc) as pool:

        for mol, xyz_confs in tqdm(pool.imap_unordered(crest, results),total=len(results),desc="Generating configurations"):

            ext_mol = Chem.Mol(mol["mol"])
            ext_mol.RemoveAllConformers()

            for conf_id, xyz_block in enumerate(xyz_confs):

                lines = xyz_block.strip().splitlines()

                n_atoms = int(lines[0])
                energy = float(lines[1])

                conf = Chem.Conformer(ext_mol.GetNumAtoms())

                for i, line in enumerate(lines[2:2 + n_atoms]):

                    _, x, y, z = line.split()

                    conf.SetAtomPosition(
                        i,
                        (float(x), float(y), float(z))
                    )

                ext_mol.AddConformer(conf, assignId=True)

                ext_mol.SetProp(
                    "ligand_id",
                    str(mol["ligand_id"])
                )

                ext_mol.SetProp(
                    "conf_id",
                    str(conf_id)
                )

                ext_mol.SetProp(
                    "energy",
                    str(energy)
                )

                writer.write(ext_mol)

    writer.close()
pass