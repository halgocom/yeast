from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import AllChem



def RdkitEtkgd(mol,ligand_id:str,etkdg_seed:int=2786245,n_conf:int=-1):
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3() #experimental torsion knowledge distance geometry (ETKDG) method 
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions =True
    params.randomSeed = etkdg_seed
    params.trackFailures = True
    params.enforceChirality = True
    params.clearConfs = True
     
    
    if n_conf == -1:
        AllChem.EmbedMolecule(mol, params)
    else:

        conf_ids = AllChem.EmbedMultipleConfs(mol,numConfs=n_conf,params=params)
        rdMolAlign.AlignMolConformers(mol)
    

    return {"id":ligand_id,"mol":mol}