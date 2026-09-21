from prody import parsePDB,saveAtoms,writePDB,AtomGroup,parseCIF
from os import getcwd,path
from rdkit import Chem
from lib import get_box

rec = ["1z95_inib.cif","1z95_r1881.cif"]
liga_name = ["inib","r1881"]
base_dir = getcwd()

for i,r in enumerate(rec):
    receptor_file = path.join(base_dir,"data",r)
    structure = parseCIF(receptor_file)
    prot = structure.select(f"chain A and not hetero")
    save_path = path.join(str(receptor_file).replace("cif","pdb"))
    writePDB(save_path, prot)

    ligand = "l01"

    ligand=structure.select(f"resname {ligand}")
    ligand_file = path.join(base_dir,"data",f"{liga_name[i]}.pdb")
    
    
    writePDB(ligand_file, ligand)
    cryst_mol = Chem.MolFromPDBFile(ligand_file,removeHs=False,sanitize=True)
    ligand_sdf_path = path.join(base_dir,"data",f"{liga_name[i]}.sdf")

    with Chem.SDWriter(ligand_sdf_path) as writer:
        writer.write(cryst_mol)

    box = get_box(cryst_mol)

    
    



    

