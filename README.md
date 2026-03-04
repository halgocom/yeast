# YeAST (Yet Another Screening Tool)

YeAST is an automated virtual screening pipeline for ligand-receptor docking, designed to streamline compound prioritization and interaction analysis.

---

## Requirements

- Install Docker CLI or Docker Desktop (Docker Desktop required for Windows)  
- Install Ubuntu WSL (Windows only)  
- Download the YeAST repository: [GitHub link](https://github.com/halgocom/yeast.git)  

---

## Installation

> **Note:** Any terminal commands on Windows must be executed inside WSL.

1. Navigate to the YeAST repository folder.  
2. Run the build script:  

```bash
./build.sh
```  

This may take a while to complete.  

3. (Docker Desktop only) Set up a volume linking the `app` folder inside the repository to `/home/screener` in the container.

---

## Usage

### Example Protocol

```python
from mk_pipeline import MeekoPipeline
from filters import LipinskiRDKit
from forcefield import RdkitEtkgd
from AutoDockScreen import calc_map, basic_docking

if __name__ == "__main__":

    clean_settings = {
        "2am9_ligand": "TES",
        "2am9_chain": "A",
        "2am9_altloc": "A",
        "pH": 6.1,
        "forcefield": "amber"
    }

    dock = MeekoPipeline(
        experiment_id="prost_canc",
        replicates=1,
        run_name="time_comp"
    )

    dock.load_targets(from_PDB=["2am9"])
    dock.target_cleanup(**clean_settings)

    dock.process_targets(
        box={"2am9": ((27.333, 1.472, 3.639), (24, 32, 16))}
    )

    # Ligands can be loaded from CSV, SDF, or compressed SDF (SDF.GZ) files
    dock.load_ligands(
        filename="cannabinoids.csv",
        sanitize=False
    )

    dock.filter_ligands(
        filter_f=LipinskiRDKit,
        allow_partial=True
    )

    dock.optimize_ligand(
        opt_f=RdkitEtkgd,
        n_conf=5,
        pH_max=6.5,
        pH_min=6.1
    )

    dock.process_ligands(output="path")

    dock.init_docking()

    dock.calculate_maps(
        map_f=calc_map,
        map_mode="o4a",
        target_id="2am9",
        target_box=True
    )

    dock.screen_single_target(
        dock_f=basic_docking,
        target_id="2am9",
        verbosity=0
    )

    dock.check_poses(target_id="2am9", redock="Testosterone")
    dock.analyze_interactions(target_id="2am9")
    dock.calc_metrics(("ki", "le"), target_id="2am9")
    dock.create_report()
```

---

### Protocol Description

The implemented protocol defines a modular and automated workflow for virtual screening on the androgen receptor ligand-binding domain (AR-LBD).  

- **`if __name__ == "__main__":`** ensures that the pipeline runs only when executed directly, preventing unintended execution during imports.  
- **`clean_settings` dictionary** allows flexible configuration of receptor parameters: ligand selection, chain, alternate locations, pH, and force field.  
- A **docking box** centered on the binding pocket ensures reproducible sampling.  
- **Ligands** are loaded from structured datasets (CSV, SDF, or SDF.GZ), filtered by Lipinski's rule of five, and optimized using the RDKit ETKDG algorithm across a defined pH range.  
- Docking maps are computed, and AutoDock-based screening is performed.  
- Post-docking tasks include pose validation, interaction analysis, and calculation of metrics like predicted inhibition constant (Ki) and ligand efficiency (LE).  
- Results can be reloaded via `from_report` for downstream analysis without repeating docking calculations.

---

### Running the Protocol

1. Place the script in `app/scripts/base`.  
2. Start the Docker container using `run.sh`.  
3. Inside the container, run:  

```bash
./screener_dock.sh
```

4. Input the protocol filename (without `.py` extension).  

This workflow allows rapid, reproducible, and automated virtual screening for ligand prioritization and binding analysis, supporting multiple ligand file formats for flexibility.
