from prody import parsePDB,saveAtoms,writePDB
import sys
from pathlib import Path

src_path = "/home/screener/files/sources/target"
pdb_in = sys.argv[1]
distance_A = float(sys.argv[2])
ligand_resname = sys.argv[3]
res_list = str(Path(src_path) / pdb_in.replace(".pdb","_pocket.txt"))
pdb_out = str(Path(src_path) / pdb_in.replace(".pdb","_pocket.pdb"))
pdb_token = str(Path(src_path) / pdb_in)





atoms_from_pdb = parsePDB(pdb_token)
prot = atoms_from_pdb.select("protein")
lig = atoms_from_pdb.select(f"resname {ligand_resname}")
lig_prot = lig + prot
writePDB(pdb_out, lig_prot)
res_selection = f"protein and same residue as within {distance_A} of resname {ligand_resname}" #select all residues within DIstance_A from ligand
residues = lig_prot.select(res_selection)
writePDB(pdb_out, residues)

header = f"Selection made with ProDy:\nLigand rename:{ligand_resname}\nSearch radius:{distance_A} A\nAtom count:{len(residues)}\n"
res_codes = list(dict.fromkeys(
    f"{a.getChid()}:{a.getResname()}{a.getResnum()}\n"
    for a in residues
))
with open(res_list,"w") as f:
    f.write(header)
    for r in res_codes:
        
        f.write(r)

