from prody import parsePDB,saveAtoms,writePDB
from urllib.request import urlretrieve
from os import getcwd,makedirs
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
#from math import round



"""
"""
def prody_fix(filepath:str):
        
    prody_string="chain A not hetero"
    
    base_pdb=parsePDB(filepath)
    file_name=filepath.split(".")[0]
    clean_pdb=f"{file_name}_clean.pdb"
    protein=base_pdb.select(prody_string)
    writePDB(clean_pdb,protein)

    return clean_pdb


def from_protein_data_bank(pdb_id:str,destination_folder:str,compressed:bool=False):

    filename = '%s.pdb' % pdb_id
    # Add .gz extension if compressed
    if compressed:
        filename = '%s.gz' % filename
    

    destination_file = f"{destination_folder}/{filename}"
    # Download the file
    url = 'https://files.rcsb.org/download/%s' % filename
    try:
        urlretrieve(url, destination_file)
    except Exception as e:
        print(e)
        
        

    if compressed:
        with gzip.open(destination_file, 'rb') as f_in:
            with open(destination_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                pass
        pass
    


    return destination_file



"""
Function to calculate vdw suface, Solvent Accessible SUrface Area, Solvent Excluded SUrface Area from Protein databank file or ID

"""



def calc_surf(protein_file:str=None,_id:str=None,clean_up:str=None,probe_radius:float=1.4,n_points:int=1000):
    _id=_id.lower()
    file_folder = f"{getcwd()}/{_id}"
    makedirs(file_folder,exist_ok=True)
    if (not protein_file) and (not _id):
        raise Exception("Missing target pdb file or target database id")

    if _id and (not protein_file):
        print(f"Downlodiang from PDB with id:{_id}")
        protein_file=from_protein_data_bank(pdb_id=_id,destination_folder=file_folder)

    print("Cleaning PDB file")
    protein=prody_fix(filepath=protein_file)


    parser = PDBParser(QUIET=1)

    struct = parser.get_structure(id=_id,file=protein)
    print("Calculating Van der Waals surface")
    VDW_sr=ShrakeRupley(probe_radius=0.01, n_points=n_points)
    VDW_sr.compute(struct, level="S")
    vdw=struct.sasa
    print(f"Calculating Solvent Accessible Surface Area\nProbe radius:{probe_radius} angstrom")
    SASA_sr=ShrakeRupley(probe_radius=probe_radius, n_points=n_points)
    SASA_sr.compute(struct, level="S")
    sasa=struct.sasa
    sesa=abs(vdw-sasa)

    vdw=vdw*100
    sasa=sasa*100
    sesa=sesa*100

    return _id,vdw,sasa,sesa
    pass






if __name__ == "__main__":
    
    _id,vdw,solv_ac,solv_ex = calc_surf(_id="bgluc_mut",protein_file="/home/screener/files/sources/target/bgluc.cif")
    with open(f"{getcwd()}/{_id}/{_id}.txt","a") as f:
        f.write(f"Surface area values for {_id} in nm2:\n vdw:{round(vdw,2)}\n SASA:{round(solv_ac,2)}\n SESA:{round(solv_ex,2)}")

    print(f"Surface area values for {_id} in nm2:\n vdw:{round(vdw,2)}\n SASA:{round(solv_ac,2)}\n SESA:{round(solv_ex,2)}")