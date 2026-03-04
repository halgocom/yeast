from typing import Iterable
import numpy as np
from tqdm import tqdm
from os import path
from urllib.request import urlretrieve
from meeko.gridbox import calc_box
from rdkit.Geometry import Point3D

def radius_of_gyration(coord_array:Iterable,padding:float):
    raise NotImplementedError
    pass


def write_pdbqt(dir:str,filename:str,pdbqt):

    pdbqt=pdbqt[0]
    with open(path.join(dir, f"{filename}.pdbqt"), "w") as f:
        f.write(pdbqt)
    return pdbqt

#Calculate normal rotation matrix
def normal_alignment(v1, v2):
    
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)
    if np.isclose(dot, -1.0):
        # rotazione 180°, scegliere asse ortogonale
        axis = np.array([1,0,0])
        if np.allclose(v1, axis):
            axis = np.array([0,1,0])
        cross = np.cross(v1, axis)
        cross /= np.linalg.norm(cross)
        K = np.array([[0,-cross[2],cross[1]],
                      [cross[2],0,-cross[0]],
                      [-cross[1],cross[0],0]])
        R = -np.eye(3) + 2*np.outer(cross,cross)
        return R
    s = np.linalg.norm(cross)
    K = np.array([[0,-cross[2],cross[1]],
                  [cross[2],0,-cross[0]],
                  [-cross[1],cross[0],0]])
    R = np.eye(3) + K + K@K*((1-dot)/(s**2))
    return R


def from_protein_data_bank(pdb_ids:Iterable[str],destination_folder:str=None,compressed:bool=False):
        paths=[]
        
        for pdb_id in tqdm(pdb_ids):
            filename = '%s.pdb' % pdb_id
            # Add .gz extension if compressed
            if compressed:
                filename = '%s.gz' % filename
            
            if destination_folder is None:
                raise RuntimeError("Missing destination folder")
            #destination_folder =  path.join(self.source_dir,"target")
            destination_file = path.join(destination_folder, filename)

            if path.exists(destination_file):
                paths.append({"target_id":pdb_id,"file":destination_file})
            else:
                print("\nDownloading receptors from ProteinDataBank\n")
                # Download the file
                url = 'https://files.rcsb.org/download/%s' % filename
                try:
                    urlretrieve(url, destination_file)
                except Exception as e:
                    print(e)
                paths.append({"id":pdb_id,"file":destination_file})

            if compressed:
                with gzip.open(destination_file, 'rb') as f_in:
                    with open(destination_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                        pass
            pass
        return paths


def get_box(mol, box_mode:str= "geometric", conf_id:int=0, padding:float=5.0): #preparation
        conf = mol.GetConformer(id=conf_id)

        coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]      #generate list with all atoms 3D coordinate
        box_dims=np.around(calc_box(coords,padding),3)
        
        #check if output is correct
        for i, values in enumerate(box_dims, start=1):
            if not hasattr(values, '__iter__'):
                raise TypeError("Returned single box value for center or dimensions")
            if len(values) != 3:
                raise ValueError("Returned less dimension or center values then expected")
            if not all(isinstance(v, float) for v in values):
                raise TypeError("Returned bad center or dimensions value types,float is expected")
        return box_dims



def coords_align(mol,align_pivot:Iterable,mol_center:Iterable,conf_id:int=-1):
    
    mol_conf=mol.GetConformer(conf_id)
    coords = np.array([list(mol_conf.GetAtomPosition(i)) for i in range(mol_conf.GetNumAtoms())])
    
    T=align_pivot-mol_center #get translation vector
    translated_coords=coords+T #apply translation vector
    for i in range(mol.GetNumAtoms()):
        x,y,z = translated_coords[i]
        mol_conf.SetAtomPosition(i,Point3D(x,y,z))
    return mol 

def map_wrapper(args,func):
    return func(*args)




def fix_pdb_columns(in_pdb, out_pdb):

    with open(in_pdb) as fin, open(out_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                serial = int(line[6:11])
                atom = line[12:16].strip()
                res = line[17:20].strip()
                chain = line[21] if line[21].strip() else "A"
                resid = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occ = float(line[54:60] or 1.00)
                temp = float(line[60:66] or 0.00)
                element = (line[76:78].strip() or atom[0]).upper()

                newline = (
                    f"{'ATOM':<6}"
                    f"{serial:>5} "
                    f"{atom:<4}"
                    f" "
                    f"{res:>3} "
                    f"{chain:1}"
                    f"{resid:>4}    "
                    f"{x:>8.3f}"
                    f"{y:>8.3f}"
                    f"{z:>8.3f}"
                    f"{occ:>6.2f}"
                    f"{temp:>6.2f}"
                    f"          "
                    f"{element:>2}\n"
                )
                fout.write(newline)
            else:
                fout.write(line)