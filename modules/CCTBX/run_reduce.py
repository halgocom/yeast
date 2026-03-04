import sys
from mmtbx.programs import reduce2
from iotbx.cli_parser import run_program



def main():
    pdb_in = sys.argv[1]
    pdb_out = sys.argv[2]


    
    phil_params=[
        "--overwrite",
        "--quiet", 
        pdb_in,
        "approach=add",
        "add_flip_movers=True",
        f"output.filename={pdb_out}"
        ]
    
    run_program(program_class=reduce2.Program,args=phil_params)
    

if __name__ == "__main__":
    main()