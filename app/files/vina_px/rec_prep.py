from prody import parsePDB,saveAtoms,writePDB,AtomGroup,parseCIF
from os import getcwd,path
from rdkit import Chem
from lib import get_box
from meeko import *
from lib import *
import json
from pathlib import Path



def prepare_target(file:str,
    target_id:str,
    altloc:str,
    box:dict,
    hydrate:bool=False,
    charge_type="gasteiger",
    config:dict=None,
    backend:str="prody",
    *args,**kwargs):

    
    outdir = str(Path(file).parent)
    outpath=path.join(outdir,f"{target_id}.pdbqt")
    box={"box_center":box[0],"box_size":box[1]}
    templates = ResidueChemTemplates.create_from_defaults() #create from defaults for now
    if config is not None:
        mk_prep = MoleculePreparation.from_config(config)
    else:
        mk_prep = MoleculePreparation(hydrate=hydrate,charge_model=charge_type)

    #load pdb with H
    if backend.lower() == "prody":
        target=parsePDB(file)

        polymer=Polymer.from_prody(
        target,
        templates,  # residue_templates, padders, ambiguous,
        mk_prep,
        #set_template,
        #delete_residues,
        allow_bad_res=True,
        #blunt_ends=True,
        #wanted_altloc=wanted_altloc,
        default_altloc=altloc
    )

    elif backend.lower() == "file":
        with open(file,"r") as pdb_file:
            pdb_string=pdb_file.read()
        
        polymer=Polymer.from_pdb_string(
        pdb_string,
        templates,  # residue_templates, padders, ambiguous,
        mk_prep,
        #set_template,
        #delete_residues,
        allow_bad_res=True,
        #blunt_ends=blunt_ends,
        #wanted_altloc=wanted_altloc,
        default_altloc=altloc,
    )
    else:
        raise ValueError(f"Invalid backend:{backend}")

    pdbqt_tuple = PDBQTWriterLegacy.write_from_polymer(polymer)
    
    pdbqt=write_pdbqt(dir=outdir,filename=f"{target_id}",pdbqt=pdbqt_tuple)

    with open(path.join(outdir,f"{target_id}.json"),"w") as json:
        json.write(polymer.to_json())

    
    return {"target_id":target_id,"path":outpath,"pdb_file":file,"json":path.join(outdir,f"{target_id}.json")} | box #return dict with target and box specs





arr=["1z95_inib","iz95_r1881"]
box=[[[-5.238, 2.604, 3.253],[19.075, 22.91, 16.975]],[[-4.889, -6.0, 0.012],[15.209, 20.878, 17.309]]]
for i,name in enumerate(arr):
    file = path.join(getcwd(),"rec_prep",f"{name}.pdb")
    preped=prepare_target(file,name,"A",box[i])

    json_file=path.join(getcwd(),"rec_prep",f"{name}.json")
    with open(json_file,"w") as f:
        json.dump(preped,f)


