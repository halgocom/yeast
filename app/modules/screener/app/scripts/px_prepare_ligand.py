

def prepare_ligand(self,mol,receptor_id:str,ligand_id:str,n_conf:int,box_mode:str="geometric",hydrate:bool=False,output:str="file",preparation_config:dict=None,etkdg_seed:int=65):

    
    print(f"Process {ligand_id} started")
    mol=Chem.AddHs(mol)
    ligand_dir=path.join(self.ligands_workpath,ligand_id)
    
    # Molecule preparation from Meeko
    if preparation_config is not None:
        preparator = MoleculePreparation.from_config(preparation_config) #load config
    else:   
        preparator = MoleculePreparation(hydrate=hydrate) #hydrate protocol

    makedirs(ligand_dir,exist_ok=True)
    
    
    
    center,size=self.get_box(mol=mol, box_mode=box_mode)
    center=np.array(center)
    size=np.array(size)



    writer = Chem.SDWriter('aspirin_confs.sdf')
    if output == "path":
        ligand_pdbqt=f"{ path.join(ligand_dir, f"{ligand_id}.pdbqt") }"
    
    # Generate conformations with 3D atom positions
    if n_conf >1:
        confs=[]
        sizes=[]
        
        #makedirs(conf_dir, exist_ok=True)

        params = AllChem.ETKDGv3() #experimental torsion knowledge distance geometry (ETKDG) method lo usiamo perchè quello di ottimizzazione stanstard di rdkit non è compatbile con mp
        params.useSmallRingTorsions = True
        params.useMacrocycleTorsions =True
        params.randomSeed = etkdg_seed
        params.trackFailures = True
        params.enforceChirality = True
        params.clearConfs = True
        conf_ids = AllChem.EmbedMultipleConfs(mol,numConfs=n_conf,params=params)
        rdMolAlign.AlignMolConformers(mol)
        for conf_id in conf_ids:
            
            Screener.coords_align(mol,center,conf_id)
            _,conf_size=self.get_box(mol=mol,box_mode=box_mode,conf_id=conf_id)
            sizes.append(conf_size)
            #centers.append(center)
            
            if hydrate: #addictional flags for hydration      
                preparator = MoleculePreparation(hydrate=hydrate) #hydrate protocol
            conf=preparator.prepare(mol,conformer_id=conf_id)
            
            conf_pdbqt = PDBQTWriterLegacy.write_string(conf[0])
            self.write_pdbqt(dir=ligand_dir,filename=f"{ligand_id}_{conf_id}",pdbqt=conf_pdbqt) #save conf pdbqt
            
            if output == "file":
                confs.append(conf_pdbqt[0])
            else:
                confs.append(f"{ path.join(ligand_dir, f"{ligand_id}_{conf_id}.pdbqt") }")
        sizes=np.array(sizes)
        opt_size=sizes.max(axis=0)
    elif n_conf == 1:
        return {"ligand_id":f"{ligand_id}","mol":mol,"base":ligand_pdbqt,"box_center":center,"box_size":size}
    else:
        raise ValueError(f"Invalid number of configurations :{n_conf}")
    
    return {"ligand_id":f"{ligand_id}","mol":mol,"base":ligand_pdbqt,"confs":tuple(confs),"box_center":center,"box_size":opt_size}
pass