from rdkit.Chem import Descriptors,MolFromSmiles

def LipinskiRDKit(smiles:str,ligand_id:str,allow_partial:bool=False):
    mol = MolFromSmiles(smiles)
    MW = Descriptors.MolWt(mol)
    HBA = Descriptors.NOCount(mol)
    HBD = Descriptors.NHOHCount(mol)
    LogP = Descriptors.MolLogP(mol)
    conditions = [MW <= 500, HBA <= 10, HBD <= 5, LogP <= 5]
    if allow_partial:
        min_count=1
    else:
        min_count=4
    
    if conditions.count(True) >= min_count:
        return  {"id":ligand_id,"mol":mol,"smiles":smiles}
    else:
        return None
    
    
"""
class LipinskiMordred(FIlter):
    pass
"""