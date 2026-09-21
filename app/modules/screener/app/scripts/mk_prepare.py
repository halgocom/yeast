from preparation import PreparationTemplate



class MeekoPreparation(PreparationTemplate):

    def process_ligands(self):
        ...
    
    def process_targets(self,box:Iterable,
        n_proc:int=None,
        hydrate:bool=False,
        preparation:Callable=mk_prepare_target,
        flexres:dict=None,
        map_mode:str="o4e",
        charge_type="geisteger",
        config:dict=None,
        ph:float=7.0,
        backend:str="prody",
        *args,
        **kwargs):


                

        #clean up and add hydrogens to receptor before preparation


        pdb_cleaner={"pdbfix":Screener.pdb_fix,"prody":Screener.prody_fix}
        PDB=pdb_cleaner.get(backend)

        if PDB is None:
            raise ValueError(f"Invalid backend\nSupported: pdbfix,prody")
        #box={}
        
        tasks=[]
        for _id in receptors.keys():
            makedirs(path.join(self.targets_workpath,_id),exist_ok=True)
            receptors[_id] = PDB(receptors[_id],self.targets_workpath,_id,hydrate)
            #add H using reduce doi.org/10.1006/jmbi.1998.2401
            #H_receptor_files = [files.replace("_clean.pdb", "_H.pdb") for files in clean_files]
        

                
            reduce=Reduce2App()
            
            results=reduce(in_files=receptors)
            print(results)
            
 
            tasks.append((self,box[_id],receptors[_id],_id,hydrate))
        
        #tasks = [(self,box[ids],files,ids,hydrate) for files,ids in zip(H_receptor_files,receptor_ids)]  #calculate task paramenters
        
        
        
        with Pool(n_proc) as pool:
            print("Start target preparation")
            results=pool.starmap(preparation, tasks)

        #sostituire questo con dataframe in shared memory
        #print(results)
        data={}
        #[{},{}] array di diz -itero array->d{},d{} se d['id'] == box['id'] 
        for column in results[0].keys():
            
            data[column] = tuple(result[column] for result in results)

        self.target_data=pd.DataFrame(data)
        self.target_data=self.target_data.set_index("id")
        print(self.target_data.head())
        
        end = time.perf_counter()
        
        print(f"Terminated in {end-start:.3f}s,avg:{(end-start)/len(results):.3f} s/receptor")

        

    def __ligand__():
        ...

    def __target__():
        pass



