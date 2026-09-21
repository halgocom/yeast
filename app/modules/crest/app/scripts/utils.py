import csv
from rdkit import Chem
import numpy as np
import multiprocessing as mp
from rdkit.Chem import rdmolfiles,AddHs,AllChem,GetPeriodicTable
from tblite import interface as tb
import subprocess
from berny import Berny, geomlib, angstrom

def optimize_xtb(ligand_id:str,smiles:str):
    
    mol = AllChem.MolFromSmiles(smiles)

    if not mol:
        return None

    base_mol = AddHs(mol)

    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    params.randomSeed = -1

    status = AllChem.EmbedMolecule(
        base_mol,
        params
    )

    xyz_block = rdmolfiles.MolToXYZBlock(base_mol)

    atomic_numbers = []
    coords = []
    pt = GetPeriodicTable()

    
    for line in xyz_block.splitlines()[2:]: #get lines and skip first 2
        sym,x,y,z = line.split()
        atomic_numbers.append(pt.GetAtomicNumber(sym))
        coords.append(list(map(float,[x,y,z])))
    
   

    atomic_numbers = np.asarray(
        atomic_numbers,
        dtype=np.int32
    )

    coords = np.asarray(
        coords,
        dtype=np.float64
    )


    geom = geomlib.Geometry(atomic_numbers, coords)

    optimizer = Berny(geom)

    geometry = next(optimizer)

    coordinates = np.asarray(
        [xyz for _, xyz in geometry],
        dtype=np.float64
    )

    calculator = tb.Calculator(
        "GFN2-xTB",
        atomic_numbers,
        coordinates * angstrom
    )

    calculator.set("verbosity", 0)

    result = calculator.singlepoint()

    energy = result["energy"]
    gradient = result["gradient"]

    optimizer.send(
        (
            energy,
            gradient / angstrom
        )
    )

    optimized_coordinates = coordinates.copy()
    previous_energy = energy

    energy_threshold = 1e-6
    gradient_threshold = 1e-4
    maxiter = 100

    for i, geometry in enumerate(optimizer, start=1):

        if i > maxiter:
            break

        coordinates = np.asarray(
            [xyz for _, xyz in geometry],
            dtype=np.float64
        )

        calculator.update(
            positions=coordinates * angstrom
        )

        result = calculator.singlepoint()

        energy = result["energy"]
        gradient = result["gradient"]

        optimized_coordinates = coordinates.copy()

        optimizer.send(
            (
                energy,
                gradient / angstrom
            )
        )

        #if delta of enenrgies between 2 consecutive iteration and gradient  are both under the respective threshold then stop optimizzaztion
        if (abs(energy - previous_energy) * 627.509474 < energy_threshold) and (np.max(np.abs(gradient)) < gradient_threshold):
            break

        previous_energy = energy

    lines = [
        str(len(atomic_numbers)),
        str(ligand_id)
    ]

    for Z, (x, y, z) in zip(
        atomic_numbers,
        optimized_coordinates
    ):
        symbol = pt.GetElementSymbol(int(Z))

        lines.append(
            f"{symbol:<2} {x: .8f} {y: .8f} {z: .8f}"
        )

    xyz_block = "\n".join(lines)

    return {
        "ligand_id": str(ligand_id),
        "xyz": xyz_block,
        "energy": energy,
        "mol":base_mol
    }

def read_csv(csv_path:str):
    with open(csv_path, newline='\n') as file:
        reader = csv.DictReader(file)
        for row in reader:
            ligand_name=row['id'].replace("-", "_").replace("/", "")

            yield (ligand_name,row["smiles"])

def read_xyz_n(path:str, n:int):
    with open(path) as f:
        for _ in range(n):
            line = f.readline()

            if not line:
                break

            n_atoms = int(line)
            comment = f.readline()
            atoms = [f.readline() for _ in range(n_atoms)]

            yield line + comment + "".join(atoms)
            
def optimize_wrapper(args):
    return optimize_xtb(*args)





def run_crest(mol:dict,crest_dir:str,algo:str,ewin:float,n_conf:int,threads:int):

    lig_dir = crest_dir / str(mol["ligand_id"])
    lig_dir.mkdir(exist_ok=True)

    xyz_path = lig_dir / f"{mol['ligand_id']}.xyz"
    log_path = lig_dir / f"{mol['ligand_id']}.log"

    with open(xyz_path, "w") as f:
        f.write(mol["xyz"])

    with open(log_path, "w") as log:
        subprocess.run(
            [
                "micromamba", "run", "-n", "crest",
                "crest", str(xyz_path),
                f"--{algo}",
                "--ewin", str(ewin),
                "-T", str(threads)
            ],
            check=True,
            cwd=lig_dir,
            stdout=log,
            stderr=subprocess.STDOUT
        )

    crest_out = lig_dir / "crest_conformers.xyz"

    xyz_confs = list(read_xyz_n(crest_out, n_conf))

    return mol, xyz_confs


def crest_wrapper(mol, crest_dir, algo, ewin, n_conf, threads):
    return run_crest(
        mol,
        crest_dir,
        algo,
        ewin,
        n_conf,
        threads
    )
