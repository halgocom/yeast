from prody import parsePDB,saveAtoms,writePDB,parseMMCIF,writeMMCIF
import sys
from pathlib import Path

src_path = "/home/screener/files/sources/target"

pdb_in = sys.argv[1]
distance_A = float(sys.argv[2])
ligand_resname = sys.argv[3]
res_list = str(Path(src_path) / pdb_in.replace(".pdb","_pocket.txt"))
pdb_out = str(Path(src_path) / pdb_in.replace(".pdb","_pocket.pdb"))
pdb_token = str(Path(src_path) / pdb_in)