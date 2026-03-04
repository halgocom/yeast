import csv 
from random import randint
import random
import gzip
import string 
from rdkit import Chem




def open_csv(filepath:str,*args):
    mols=set()
    with open(filepath, newline='\n') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row['id'].isalnum(): #da cambiare
                ligand_name=row['id']
            else:
                ligand_name=''.join(random.SystemRandom().choice(string.ascii_uppercase) for _ in range(6))
                #ids.add(ligand_name) #generate random alphanumerical id
            mol=Chem.MolFromSmiles(SMILES=row['smiles'])
            if mol is None:
                continue
            
            #print(ligand_name)
            mol.SetProp("id",ligand_name)
            mols.add(mol)
    pass
    
    return mols


def sdf_from_supplier(suppl,basepath:str):
        mols=set()
        for mol in supplier:
            if mol is None:
                continue
            
            cid=''.join(random.choice(chars) for _ in range(size)) #assign random id
            mol.SetProp("id",cid)
            ligand_path=path.join(base_path,f"{cid}")
            makedirs(ligand_path,exist_ok=True)
            outsdf_path=f"{ligand_path}/{cid}.sdf"
            mols.add(mol)
            writer = Chem.SDWriter(outsdf_path)
            writer.write(mol)
            writer.close()
        return mols


def open_sdf(filepath:str,filetype:str,base_path:str):

    if filetype == "sdf":
        file=open(filepath,"rb")
    else:
        file=gzip.open(filepath, 'rb')
    
    with Chem.ForwardSDMolSupplier(file) as supply:
        mols,ids=sdf_from_supplier(suppl=supply,basepath=base_path)    
    
    return mols,ids

def fix_ligands(smiles: str, _id: str):
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return None

        mol.SetProp("id", _id)

        smiles =Chem.MolToSmiles(
            mol,
            canonical=True,
            isomericSmiles=True
        )
        
        return {"id":_id,"smiles":smiles}
    except Exception:
        
        return None