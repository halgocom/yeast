from typing import Iterable
import numpy as np
from tqdm import tqdm
from os import path
from urllib.request import urlretrieve
from meeko.gridbox import calc_box
from rdkit.Geometry import Point3D
from pathlib import Path

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


def from_protein_data_bank(pdb_id:str,destination_folder:str=None):
    paths=[]
    
    pdb_id = pdb_id.split(".")[0]


    filename = '%s.pdb' % pdb_id

    if destination_folder is None:
        raise RuntimeError("Missing destination folder")

    destination_file = str(Path(destination_folder) /filename)

    # Download the file
    url = 'https://files.rcsb.org/download/%s' % filename

    urlretrieve(url, str(destination_file))

    return {"target_id":pdb_id,"file":destination_file}


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

def box_to_pdb_string(box_center, npts, spacing=0.375):
    eol="\n"
    step_x = int(npts[0] / 2.0) * spacing
    step_y = int(npts[1] / 2.0) * spacing
    step_z = int(npts[2] / 2.0) * spacing
    center_x, center_y, center_z = box_center
    corners = []
    corners.append([center_x - step_x, center_y - step_y, center_z - step_z] ) # 1
    corners.append([center_x + step_x, center_y - step_y, center_z - step_z] ) # 2
    corners.append([center_x + step_x, center_y + step_y, center_z - step_z] ) # 3
    corners.append([center_x - step_x, center_y + step_y, center_z - step_z] ) # 4
    corners.append([center_x - step_x, center_y - step_y, center_z + step_z] ) # 5
    corners.append([center_x + step_x, center_y - step_y, center_z + step_z] ) # 6
    corners.append([center_x + step_x, center_y + step_y, center_z + step_z] ) # 7
    corners.append([center_x - step_x, center_y + step_y, center_z + step_z] ) # 8

    count = 1
    res = "BOX"
    chain = "X"
    pdb_out = ""
    line = "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00          %1s" + eol
    for idx in range(len(corners)):
        x = corners[idx][0]
        y = corners[idx][1]
        z = corners[idx][2]
        pdb_out += line % (count, "Ne", res, chain, idx+1, x, y, z, "Ne")
        count += 1

    # center
    pdb_out += line % (count+1, "Xe", res, chain, idx+1, center_x, center_y, center_z, "Xe")

    pdb_out += "CONECT    1    2" + eol
    pdb_out += "CONECT    1    4" + eol
    pdb_out += "CONECT    1    5" + eol
    pdb_out += "CONECT    2    3" + eol
    pdb_out += "CONECT    2    6" + eol
    pdb_out += "CONECT    3    4" + eol
    pdb_out += "CONECT    3    7" + eol
    pdb_out += "CONECT    4    8" + eol
    pdb_out += "CONECT    5    6" + eol
    pdb_out += "CONECT    5    8" + eol
    pdb_out += "CONECT    6    7" + eol
    pdb_out += "CONECT    7    8" + eol
    return pdb_out



