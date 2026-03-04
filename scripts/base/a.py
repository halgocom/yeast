"""
from rdkit import Chem
from rdkit.Chem import AllChem
"""
smiles ="CCCCCC1=CC(=C2[C@@H]3C=C(CC[C@H]3C(OC2=C1)(C)C)C)O"
from rdkit import Chem
from rdkit.Chem import EnumerateStereoisomers
molecule = Chem.MolFromSmiles(smiles)
options = EnumerateStereoisomers.StereoEnumerationOptions(unique=True, tryEmbedding=True)
isomers = tuple(EnumerateStereoisomers.EnumerateStereoisomers(
      molecule,
      options=options)
)
for smiles in isomers:
  print(Chem.MolToSmiles(smiles, isomericSmiles=True))
"""
import glob
import os
import numpy as np

from numpy import std
from scipy import stats
import matplotlib as mlp
from typing import Iterable
import pandas as pd
import csv
import math as m 

# Directory in cui si trova lo script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Cerca i file .sdf in quella directory
sdf_files = glob.glob(os.path.join(script_dir, "*.sdf"))
print(f"Trovati {len(sdf_files)} file SDF")

writer = Chem.SDWriter(os.path.join(script_dir, "merged_output.sdf"))
count = 0

for sdf_file in sdf_files:
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        if mol is not None:
            writer.write(mol)
            count += 1

writer.close()
print(f"Mergiate {count} molecole.")



"""
























































































































def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def luminescence_norm(blank_lum:float,data:Iterable=None,concentrations:Iterable=None):
    if (data is None) or (concentrations is None):
        raise Exception("Missing experimental data")
    log10_conc = [m.log10(conc) for conc in concentrations]
    norm_lum = [lum/blank_lum for lum in data]
    return log10_conc,norm_data


def EC50_lum_varioskan(ligand_id:str,receptor_id:str,concentrations:Iterable,outlier:dict=None,filename:str=None,conc_unit:str="%w/v",plot:bool=True,replicates=3,rep_separator:str="_",with_fold_induction:bool=True):
    if filename is None:
        raise Exception("Missing experimental data filename in screener/data")
    
    lum_data=pd.read_csv(filename,dtype=np.float64)
    print(lum_data.head(16))
    data_dict={}
    mean=[]
    dev_st=[]
    dev_st_perc=[]
    dev_avg=[]
    norm_error=[]
    norm_rlu=[]
    norm_error_blank=[]

    blank_max=[lum_data[f"Blank_{i}"].max() for i in range(1,replicates+1)]
    mean.append(sum(blank_max)/len(blank_max))
    dev_st.append(std(blank_max))
    dev_st_perc.append((dev_st[0]/mean[0])*100)
    dev_avg.append(dev_st[0]/len(blank_max))
    norm_rlu.append((mean[0]/mean[0])*100)
    norm_error.append(m.sqrt(m.pow(dev_st_perc[0]/100,2)+m.pow(dev_st_perc[0]/100,2))*mean[0])
    #Mettere controllo per none su outliers
    for i,conc in enumerate(concentrations):
        max_rlu=[]
        for j in range(1,replicates+1):
            conc_outliers=outlier.get(str(conc))

            try:
                if (outlier is not None) and (j in conc_outliers): #skips replicate denoted as outlier
                    continue
            except Exception:
                pass 
            if not with_fold_induction:
                lum_data[f"{conc}_{j}"] = lum_data[f"{conc}_{j}"]-lum_data[f"Blank_{j}"] #blank subtraction
            rlu=lum_data[f"{conc}_{j}"].max() #find max rlu
            max_rlu.append(rlu)
        
        mean.append(sum(max_rlu)/len(max_rlu)) #calculate max rlu value
        dev_st.append(std(max_rlu)) #calculete dev st
        dev_st_perc.append((dev_st[i+1]/mean[0])*100) #calculate dev st %
        dev_avg.append(dev_st[i+1]/len(max_rlu)) #calculate dev st avg
    norm_rlu=[(avg/max(mean)*100) for avg in mean] #normazation to 100 of rlu data
    
    error_ctr=m.pow(dev_st_perc[0]/100,2)
    norm_error=[m.sqrt(error_ctr+m.pow(dev_st_perc[i]/100,2))*norm_rlu[i] for i in range(len(norm_rlu))] #calculate error on normalized RLU
        
    concentrations.insert(0,0)
    if with_fold_induction:
        fi=[avg/mean[0] for avg in mean]
        fi_error=[m.sqrt(error_ctr+m.pow(dev_st_perc[i]/100,2))*fi for i,fi in enumerate(fi)]
        data = {f"Concentration ({conc_unit})":concentrations,"Average RLU":mean,"Normalized RLU":norm_rlu,"Standard Deviation":dev_st,"Relative Standard Deviation %":dev_st_perc,"Normalized Error":norm_error,"Fold Induction":fi,"Fold Induction Error":fi_error}
    else:
        data = {f"Concentration ({conc_unit})":concentrations,"Average RLU":mean,"Normalized RLU":norm_rlu,"Standard Deviation":dev_st,"Relative Standard Deviation %":dev_st_perc,"Normalized Error":norm_error}
   
    data_df = pd.DataFrame(data)
    print(data_df.head(8))
    data_df.to_csv("/home/screener/data/tnfhar.csv",sep=",",encoding="utf-8",index=False)
            


    #print(lum_data.head())
    
   
    

    
    




EC50_lum_varioskan(ligand_id="TNFalpha",receptor_id="hAR",filename="/home/screener/data/TNF_hAR.csv",conc_unit="ng/mL",concentrations=[0.1,0.5,1,2.5,5,10,20],outlier={"20":[3]})
