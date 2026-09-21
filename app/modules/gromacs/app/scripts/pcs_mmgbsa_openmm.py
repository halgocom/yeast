"""
Single-trajectory MD + MM/GBSA workflow
=======================================

Workflow:
    docked SDF
      -> ligand PDB
      -> ACPYPE/GAFF2 topology
      -> protein + retained crystallographic waters
      -> OpenMM solvation/ions
      -> OpenMM explicit-solvent MD
      -> MDTraj trajectory analysis
      -> OpenMM OBC2/GBSA single-point energies
      -> MM/GBSA statistics and plots

Important:
- The MD is performed with OpenMM.
- MDTraj is used for trajectory analysis.
- ACPYPE is used for ligand GAFF2/AM1-BCC-style GROMACS parameter generation.
- ParmEd is used as the bridge between the GROMACS/ACPYPE topology and OpenMM,
  and to construct the explicit/implicit-solvent OpenMM systems.
- MM/GBSA is single-trajectory: receptor, ligand and complex coordinates all
  come from the same complex trajectory.
- Entropy (-TΔS) is NOT calculated here. The reported quantity is therefore
  ΔG_MMGBSA = ΔE_MM + ΔG_GB + ΔG_SA, not a full ΔG_bind including entropy.

The script deliberately keeps all settings in the configuration section.
"""

from pathlib import Path
import glob
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mdtraj as md
import parmed as pmd
from parmed import unit as pmd_unit

from rdkit import Chem
from prody import parsePDB, writePDB

from acpype.topol import ACTopol

from openmm import LangevinMiddleIntegrator,MonteCarloBarostat,Platform,XmlSerializer
from openmm import unit
from openmm.app import PDBFile,ForceField,Modeller,Simulation,PME,NoCutoff,HBonds,OBC2,DCDReporter,StateDataReporter


# CONFIGURATION


EXPERIMENT_NAME = "prost_canc"
RUN_NAME = "vina_new"

MD_PATH = Path(
    f"/home/gromacs/shared_files/{EXPERIMENT_NAME}/{RUN_NAME}/md"
)

LIGAND_ID = "5FAKB48"
TARGET_ID = "2am9"

TARGET_PATH = Path(
    f"/home/gromacs/shared_files/{EXPERIMENT_NAME}/{RUN_NAME}/"
    f"targets/{TARGET_ID}/{TARGET_ID}_fixed.pdb"
)

DOCKING_PATH = Path(
    f"/home/gromacs/shared_files/{EXPERIMENT_NAME}/{RUN_NAME}/results"
)

LIGAND_PATH = Path(
    f"{DOCKING_PATH}/{TARGET_ID}/replicate 1/{LIGAND_ID}/"
    f"{TARGET_ID}_{LIGAND_ID}_1_replicate1.sdf"
)

ORIG_TARGET = Path(
    f"/home/gromacs/shared_files/sources/target/{TARGET_ID}.pdb"
)

LIGAND_OUTDIR = MD_PATH / LIGAND_ID
RESULTS_DIR = LIGAND_OUTDIR / "mmgbsa"
PLOT_DIR = RESULTS_DIR / "plots"


# Docking-box water selection.
POCKET_CENTER = [1,2,3]

POCKET_CENTER = np.array(POCKET_CENTER, dtype=float)

POCKET_SIZE = np.array(POCKET_SIZE, dtype=float)


# Force fields.


PROTEIN_FF = "amber99sb-ildn.xml"
WATER_FF = "tip3p.xml"

# ACPYPE/Antechamber ligand parameters.
LIGAND_ATOM_TYPE = "gaff2"
LIGAND_CHARGE_TYPE = "bcc"
LIGAND_CHARGE = 0

# Physiological ionic strength used during explicit-solvent preparation and
# as the ionic strength for the implicit GB calculation.
IONIC_STRENGTH = 0.15 * unit.molar

TEMPERATURE = 300 * unit.kelvin
PRESSURE = 1 * unit.atmosphere

# 2 fs is deliberately used rather than 4 fs because no hydrogen-mass
# repartitioning is performed here.
TIMESTEP = 0.002 * unit.picoseconds
FRICTION = 1.0 / unit.picoseconds


# MD protocol.


NVT_EQUIL_STEPS = 50_000       # 100 ps
NPT_EQUIL_STEPS = 250_000      # 500 ps
PRODUCTION_STEPS = 5_000_000   # 10 ns

TRAJ_INTERVAL = 5_000          # 10 ps
ENERGY_INTERVAL = 5_000        # 10 ps

# Analyse only the last fraction of the production trajectory.
ANALYSIS_FRACTION = 0.50

# To avoid unnecessarily expensive single-point calculations, analyse every
# Nth saved production frame.
MMGBSA_FRAME_STRIDE = 1


OPENMM_PLATFORM = "CPU"


GB_SOLUTE_DIELECTRIC = 1.0
GB_SOLVENT_DIELECTRIC = 80.0


# BASIC VALIDATION


MD_PATH.mkdir(parents=True, exist_ok=True)
LIGAND_OUTDIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

if POCKET_CENTER.size != 3 or POCKET_SIZE.size != 3:
    raise ValueError(
        "POCKET_CENTER and POCKET_SIZE must each contain exactly "
        "three coordinates in Å."
    )

if np.any(POCKET_SIZE <= 0):
    raise ValueError("POCKET_SIZE must contain positive full box dimensions.")

if not TARGET_PATH.exists():
    raise FileNotFoundError(TARGET_PATH)

if not LIGAND_PATH.exists():
    raise FileNotFoundError(LIGAND_PATH)

if not ORIG_TARGET.exists():
    raise FileNotFoundError(ORIG_TARGET)

# FILE NAMES


LIGAND_PDB = MD_PATH / f"{LIGAND_ID}.pdb"
POCKET_WATER_PDB = MD_PATH / f"{TARGET_ID}_W.pdb"
COMPLEX_DRY_PDB = MD_PATH / f"{TARGET_ID}_{LIGAND_ID}.pdb"
COMPLEX_SOLVATED_PDB = MD_PATH / f"{TARGET_ID}_{LIGAND_ID}_solvated.pdb"

LIGAND_ACPYPE_DIR = LIGAND_OUTDIR / f"{LIGAND_ID}.acpype"

# ACPYPE normally writes these names for GROMACS output.
LIGAND_TOP = LIGAND_ACPYPE_DIR / f"{LIGAND_ID}_GMX.top"
LIGAND_GRO = LIGAND_ACPYPE_DIR / f"{LIGAND_ID}_GMX.gro"

TRAJECTORY_DCD = RESULTS_DIR / "production.dcd"
TRAJECTORY_PDB = RESULTS_DIR / "trajectory_topology.pdb"
ENERGY_CSV = RESULTS_DIR / "mmgbsa_per_frame.csv"
SUMMARY_CSV = RESULTS_DIR / "mmgbsa_summary.csv"


# HELPER FUNCTIONS


def select_openmm_platform(name):
    """Return an OpenMM Platform by name."""
    return Platform.getPlatformByName(name)


def write_rdkit_ligand_pdb(mol, output_path, residue_name):

    if mol is None:
        raise ValueError("The SDF did not contain a valid molecule.")

    if mol.GetNumConformers() == 0:
        raise ValueError("The ligand has no 3D conformer.")

    for atom in mol.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol():>2}")
        info.SetResidueName(residue_name)
        info.SetResidueNumber(1)
        info.SetChainId("L")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)

    Chem.MolToPDBFile(mol, str(output_path))


def create_ligand_topology():
    """
    Generate a GAFF2 ligand topology through ACPYPE's Python interface.
    """
    ligand_mol = Chem.SDMolSupplier(
        str(LIGAND_PATH),
        removeHs=False,
    )[0]

    write_rdkit_ligand_pdb(
        ligand_mol,
        LIGAND_PDB,
        LIGAND_ID[:3].upper(),
    )

    os.chdir(str(LIGAND_OUTDIR))


    # Wang et al., J. Comput. Chem. 2004, 25, 1157-1174.
    # GAFF was developed to provide parameters for drug-like organic molecules while remaining compatible with AMBER biomolecular FFs.

    gmx_ligand = ACTopol(
        verbose=False,
        inputFile=str(LIGAND_PDB),
        chargeType=LIGAND_CHARGE_TYPE,
        chargeVal=LIGAND_CHARGE,
        atomType=LIGAND_ATOM_TYPE,
        basename=LIGAND_ID,
        chiral=True,
        outTopol="gmx",
    )

    gmx_ligand.createACTopol()
    gmx_ligand.createMolTopol()
    gmx_ligand.molTopol.writeGromacsTopolFiles()


    candidates = [
        LIGAND_TOP,
        LIGAND_OUTDIR / f"{LIGAND_ID}_GMX.top",
        LIGAND_OUTDIR / f"{LIGAND_ID}.acpype" / f"{LIGAND_ID}_GMX.top",
    ]

    if not any(p.exists() for p in candidates):
        found = list(LIGAND_OUTDIR.rglob("*_GMX.top"))
        if found:
            return Path(found[0])

        raise FileNotFoundError(
            "ACPYPE finished but no *_GMX.top file was found under "
            f"{LIGAND_OUTDIR}"
        )

    return next(p for p in candidates if p.exists())


def prepare_pocket_waters():

    original = parsePDB(str(ORIG_TARGET))

    high = POCKET_CENTER + POCKET_SIZE / 2.0
    low = POCKET_CENTER - POCKET_SIZE / 2.0

    selection = (
        f"water and "
        f"{low[0]} <= x <= {high[0]} and "
        f"{low[1]} <= y <= {high[1]} and "
        f"{low[2]} <= z <= {high[2]}"
    )

    waters = original.select(selection)

    if waters is None or waters.numAtoms() == 0:
        print("WARNING: no crystallographic waters were found in the docking box.")
        return None

    writePDB(str(POCKET_WATER_PDB), waters)
    return POCKET_WATER_PDB


def build_dry_complex(pocket_water_path):
    """
    Build target + docked ligand + selected crystallographic waters.
    """
    target = parsePDB(str(TARGET_PATH))
    ligand = parsePDB(str(LIGAND_PDB))

    complex_atoms = target + ligand

    if pocket_water_path is not None:
        waters = parsePDB(str(pocket_water_path))
        complex_atoms = complex_atoms + waters

    writePDB(str(COMPLEX_DRY_PDB), complex_atoms)
    return COMPLEX_DRY_PDB


def solvate_with_openmm(dry_complex_path):
    """
    Add missing protein hydrogens and an explicit TIP3P/0.15 M solvent box.

    The ligand is present during solvation, but it is NOT parameterized by the
    protein force field. Parameterization is done separately with ACPYPE.
    """
    pdb = PDBFile(str(dry_complex_path))

    # AMBER99SB-ILDN + TIP3P is used here as the biomolecular/solvent model.

    forcefield = ForceField(PROTEIN_FF, WATER_FF)
    # First process the protein/water part alone so that the protein FF does not need to recognize the non-standard ligand residue.
    ligand_resname = LIGAND_ID[:3].upper()

    ligand_atoms = [
        atom
        for atom in pdb.topology.atoms()
        if atom.residue.name.strip().upper() == ligand_resname
    ]

    if not ligand_atoms:
        raise RuntimeError(
            f"Could not find ligand residue '{ligand_resname}' in "
            f"{dry_complex_path}"
        )

    base_modeller = Modeller(pdb.topology, pdb.positions)
    base_modeller.delete(ligand_atoms)


    base_modeller.addHydrogens(forcefield, pH=7.4)

    # Add the ligand back with its original docked coordinates.
    ligand_pdb = PDBFile(str(LIGAND_PDB))
    base_modeller.add(ligand_pdb.topology, ligand_pdb.positions)


    # Explicit TIP3P solvent and neutralizing ions are added, plus NaCl toreach the requested ionic strength.

    base_modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=1.0 * unit.nanometer,
        ionicStrength=IONIC_STRENGTH,
        neutralize=True,
        positiveIon="Na+",
        negativeIon="Cl-",
    )

    with open(COMPLEX_SOLVATED_PDB, "w") as handle:
        PDBFile.writeFile(
            base_modeller.topology,
            base_modeller.positions,
            handle,
            keepIds=True,
        )

    return base_modeller, forcefield


def parameterize_explicit_system(modeller, forcefield, ligand_top):

    ligand_resname = LIGAND_ID[:3].upper()

    ligand_atoms = [
        atom
        for atom in modeller.topology.atoms()
        if atom.residue.name.strip().upper() == ligand_resname
    ]

    if not ligand_atoms:
        raise RuntimeError("Ligand atoms disappeared during solvation.")

    # Make a copy and remove the ligand so amber99sb-ildn only sees residues
    # for which it has templates.
    receptor_modeller = Modeller(
        modeller.topology,
        modeller.positions,
    )
    receptor_modeller.delete(ligand_atoms)


    receptor_system = forcefield.createSystem(
        receptor_modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=5e-5,
    )

    # Convert the OpenMM topology + System into a fully parametrized ParmEd
    # Structure. ParmEd can then be combined with the ACPYPE/GROMACS ligand.
    receptor_parm = pmd.openmm.load_topology(
        receptor_modeller.topology,
        system=receptor_system,
        xyz=np.asarray(
            receptor_modeller.positions.value_in_unit(unit.angstrom)
        ),
    )

    ligand_parm = pmd.load_file(
        str(ligand_top),
        xyz=str(LIGAND_GRO),
    )

    # The ligand topology must contain exactly the same number of atoms as
    # the docked ligand PDB.
    if len(ligand_parm.atoms) != len(ligand_atoms):
        raise ValueError(
            "ACPYPE ligand topology atom count does not match the docked "
            "ligand coordinates: "
            f"{len(ligand_parm.atoms)} vs {len(ligand_atoms)}"
        )

    complex_parm = receptor_parm + ligand_parm

    # Preserve the periodic box generated by OpenMM.
    complex_parm.box = receptor_parm.box

    # Save a self-contained GROMACS representation for provenance/reuse.
    complex_parm.save(
        str(RESULTS_DIR / "complex.gro"),
        overwrite=True,
    )
    complex_parm.save(
        str(RESULTS_DIR / "complex.top"),
        overwrite=True,
    )

    return complex_parm


def choose_platform():
    return select_openmm_platform(OPENMM_PLATFORM)


def make_explicit_md_system(complex_parm):

    # Essmann et al., J. Chem. Phys. 1995, 103, 8577-8593.
    # PME is used for long-range electrostatics in the periodic explicit-
    # solvent MD system.

    system = complex_parm.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * pmd_unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=5e-5,
        removeCMMotion=True,
    )

    return system


def run_md(complex_parm):

    system = make_explicit_md_system(complex_parm)

    integrator = LangevinMiddleIntegrator(
        TEMPERATURE,
        FRICTION,
        TIMESTEP,
    )

    platform = choose_platform()

    simulation = Simulation(
        complex_parm.topology,
        system,
        integrator,
        platform,
    )

    simulation.context.setPositions(complex_parm.positions)

    # Zhang et al., J. Phys. Chem. A 2019, 123, 6056-6079,
    # DOI 10.1021/acs.jpca.9b02771.


    print("Minimization...")
    simulation.minimizeEnergy(maxIterations=10_000)

    print("NVT equilibration...")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)

    simulation.reporters.append(
        StateDataReporter(
            str(RESULTS_DIR / "nvt.log"),
            ENERGY_INTERVAL,
            step=True,
            time=True,
            temperature=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
        )
    )

    simulation.step(NVT_EQUIL_STEPS)

   --------
    # MonteCarloBarostat is used only for NPT equilibration
    # Pressure coupling should not be confused with the GB calculation the trajectory is generated in explicit solvent, and GBSA is applied afterward to stripped solute snapshots.
   --------
    system.addForce(
        MonteCarloBarostat(
            PRESSURE,
            TEMPERATURE,
            25,
        )
    )

    print("NPT equilibration...")
    simulation.context.reinitialize(preserveState=True)

    simulation.reporters.clear()
    simulation.reporters.append(
        StateDataReporter(
            str(RESULTS_DIR / "npt.log"),
            ENERGY_INTERVAL,
            step=True,
            time=True,
            temperature=True,
            potentialEnergy=True,
            volume=True,
        )
    )

    simulation.step(NPT_EQUIL_STEPS)

    print("Production MD...")
    simulation.reporters.clear()

    simulation.reporters.append(
        DCDReporter(
            str(TRAJECTORY_DCD),
            TRAJ_INTERVAL,
            enforcePeriodicBox=True,
        )
    )

    simulation.reporters.append(
        StateDataReporter(
            str(RESULTS_DIR / "production.log"),
            ENERGY_INTERVAL,
            step=True,
            time=True,
            temperature=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            volume=True,
        )
    )

    simulation.step(PRODUCTION_STEPS)

    # Save final coordinates.
    final_state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
    )

    with open(RESULTS_DIR / "final.pdb", "w") as handle:
        PDBFile.writeFile(
            simulation.topology,
            final_state.getPositions(),
            handle,
            keepIds=True,
        )

    return simulation


def build_analysis_topology(complex_parm):
    """
    Write a PDB containing exactly the atom order used in the trajectory.
    MDTraj uses this PDB only as topology information.
    """
    complex_parm.save(
        str(TRAJECTORY_PDB),
        overwrite=True,
    )
    return TRAJECTORY_PDB


def load_production_trajectory(topology_path):
    """
    Load the DCD with MDTraj and select the production analysis window.
    """
    traj = md.load(
        str(TRAJECTORY_DCD),
        top=str(topology_path),
    )

    if traj.n_frames < 10:
        raise RuntimeError(
            f"Only {traj.n_frames} frames were found in the production DCD."
        )

    first = int(traj.n_frames * (1.0 - ANALYSIS_FRACTION))
    analysis = traj[first::MMGBSA_FRAME_STRIDE]

    return traj, analysis



# McGibbon et al., Biophys. J. 2015, 109, 1528-1532,
# DOI 10.1016/j.bpj.2015.08.015.
# MDTraj is used here for RMSD/RMSF/Rg/SASA trajectory analysis.

def identify_indices(traj):
    """
    Identify protein, ligand and ligand-heavy-atom indices from the MDTraj
    topology.
    """
    ligand_resname = LIGAND_ID[:3].upper()

    ligand_atoms = [
        atom.index
        for atom in traj.topology.atoms
        if atom.residue.name.strip().upper() == ligand_resname
    ]

    if not ligand_atoms:
        raise RuntimeError(
            f"MDTraj could not identify ligand residue '{ligand_resname}'."
        )

    protein_atoms = traj.topology.select("protein")

    ligand_heavy = np.array([
        i for i in ligand_atoms
        if traj.topology.atom(i).element is not None
        and traj.topology.atom(i).element.atomic_number != 1
    ], dtype=int)

    protein_backbone = traj.topology.select(
        "protein and backbone"
    )

    protein_ca = traj.topology.select(
        "protein and name CA"
    )

    return (
        np.asarray(protein_atoms, dtype=int),
        np.asarray(ligand_atoms, dtype=int),
        ligand_heavy,
        np.asarray(protein_backbone, dtype=int),
        np.asarray(protein_ca, dtype=int),
    )


def plot_trajectory_analysis(traj, analysis, indices):
    """
    Generate trajectory-quality-control plots with MDTraj.
    """
    (
        protein_atoms,
        ligand_atoms,
        ligand_heavy,
        protein_backbone,
        protein_ca,
    ) = indices

    time_ns = analysis.time / 1000.0


    # Protein backbone RMSD after C-alpha alignment.

    aligned = analysis.superpose(
        analysis,
        0,
        atom_indices=protein_ca,
    )

    protein_rmsd = md.rmsd(
        aligned,
        aligned[0],
        frame=0,
        atom_indices=protein_backbone,
    ) * 10.0

    plt.figure(figsize=(8, 5))
    plt.plot(time_ns, protein_rmsd)
    plt.xlabel("Time (ns)")
    plt.ylabel("Protein backbone RMSD (Å)")
    plt.title("Protein backbone RMSD")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "protein_rmsd.png", dpi=300)
    plt.close()


    # Ligand heavy-atom RMSD in the protein-aligned trajectory.

    ligand_rmsd = md.rmsd(
        aligned,
        aligned[0],
        frame=0,
        atom_indices=ligand_heavy,
    ) * 10.0

    plt.figure(figsize=(8, 5))
    plt.plot(time_ns, ligand_rmsd)
    plt.xlabel("Time (ns)")
    plt.ylabel("Ligand heavy-atom RMSD (Å)")
    plt.title("Ligand RMSD after protein alignment")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "ligand_rmsd.png", dpi=300)
    plt.close()


    # Protein C-alpha RMSF.

    rmsf = md.rmsf(
        aligned,
        aligned[0],
        atom_indices=protein_ca,
    ) * 10.0

    ca_residues = [
        aligned.topology.atom(i).residue.index
        for i in protein_ca
    ]

    plt.figure(figsize=(10, 5))
    plt.plot(ca_residues, rmsf)
    plt.xlabel("Residue index")
    plt.ylabel("Cα RMSF (Å)")
    plt.title("Protein Cα RMSF")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "protein_rmsf.png", dpi=300)
    plt.close()


    # Ligand radius of gyration.

    ligand_rg = md.compute_rg(
        analysis.atom_slice(ligand_heavy)
    ) * 10.0

    plt.figure(figsize=(8, 5))
    plt.plot(time_ns, ligand_rg)
    plt.xlabel("Time (ns)")
    plt.ylabel("Ligand radius of gyration (Å)")
    plt.title("Ligand compactness")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "ligand_rg.png", dpi=300)
    plt.close()


    # Solvent-accessible surface area.

    ligand_sasa = md.shrake_rupley(
        analysis.atom_slice(ligand_heavy),
        mode="atom",
        probe_radius=0.14,
    ).sum(axis=1) * 100.0

    plt.figure(figsize=(8, 5))
    plt.plot(time_ns, ligand_sasa)
    plt.xlabel("Time (ns)")
    plt.ylabel("Ligand SASA (Å²)")
    plt.title("Ligand solvent-accessible surface area")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "ligand_sasa.png", dpi=300)
    plt.close()

    return {
        "protein_rmsd_A": protein_rmsd,
        "ligand_rmsd_A": ligand_rmsd,
        "protein_rmsf_A": rmsf,
        "ligand_rg_A": ligand_rg,
        "ligand_sasa_A2": ligand_sasa,
    }



# Gohlke & Case / Miller et al. end-state MM-PB(GB)SA framework:
# ΔG_bind = G_complex - G_receptor - G_ligand, evaluated from snapshots of
# an MD ensemble. The single-trajectory protocol uses the same complex
# trajectory for all three components.

# Miller et al., J. Chem. Theory Comput. 2012, 8, 3314-3321,
# DOI 10.1021/ct300418h.

def create_vacuum_system(structure):
    """
    Create a gas-phase MM system.

    All covalent terms remain present. No solvent is included.
    """
    return structure.createSystem(
        nonbondedMethod=NoCutoff,
        constraints=None,
        rigidWater=False,
        implicitSolvent=None,
        removeCMMotion=False,
    )


# --------------------------------------------------------------------------
# Onufriev, Bashford & Case, Proteins 2004, 55, 383-394.
# OBC2 is the AMBER igb=5 GBSA model; OpenMM exposes it directly as OBC2.
# --------------------------------------------------------------------------
def create_gbsa_system(structure):
    """
    Create the implicit-solvent OBC2/GBSA system used for the GBSA term.

    OBC2 corresponds to AMBER igb=5.
    """
    return structure.createSystem(
        nonbondedMethod=NoCutoff,
        constraints=None,
        rigidWater=False,
        implicitSolvent=OBC2,
        implicitSolventKappa=None,
        implicitSolventSaltConc=IONIC_STRENGTH,
        soluteDielectric=GB_SOLUTE_DIELECTRIC,
        solventDielectric=GB_SOLVENT_DIELECTRIC,
        useSASA=True,
        removeCMMotion=False,
    )


def clone_system(system):
    return XmlSerializer.deserialize(
        XmlSerializer.serialize(system)
    )


def modify_nonbonded(system, mode):
    """
    Return a system in which only the requested nonbonded contribution is
    retained, while bonded terms remain untouched.

    mode:
        "bonded"       -> all charges and LJ epsilons set to zero
        "electrostatic"-> LJ epsilons set to zero
        "vdw"          -> charges set to zero
    """
    new_system = clone_system(system)

    nonbonded = None

    for force in new_system.getForces():
        if force.__class__.__name__ == "NonbondedForce":
            nonbonded = force
            break

    if nonbonded is None:
        raise RuntimeError(
            "No OpenMM NonbondedForce was found in the vacuum system."
        )

    for i in range(nonbonded.getNumParticles()):
        q, sigma, epsilon = nonbonded.getParticleParameters(i)

        if mode == "bonded":
            nonbonded.setParticleParameters(
                i,
                0 * q,
                sigma,
                0 * epsilon,
            )
        elif mode == "electrostatic":
            nonbonded.setParticleParameters(
                i,
                q,
                sigma,
                0 * epsilon,
            )
        elif mode == "vdw":
            nonbonded.setParticleParameters(
                i,
                0 * q,
                sigma,
                epsilon,
            )
        else:
            raise ValueError(mode)

    for i in range(nonbonded.getNumExceptions()):
        p1, p2, chargeprod, sigma, epsilon = (
            nonbonded.getExceptionParameters(i)
        )

        if mode == "bonded":
            nonbonded.setExceptionParameters(
                i,
                p1,
                p2,
                0 * chargeprod,
                sigma,
                0 * epsilon,
            )
        elif mode == "electrostatic":
            nonbonded.setExceptionParameters(
                i,
                p1,
                p2,
                chargeprod,
                sigma,
                0 * epsilon,
            )
        elif mode == "vdw":
            nonbonded.setExceptionParameters(
                i,
                p1,
                p2,
                0 * chargeprod,
                sigma,
                epsilon,
            )

    return new_system


def make_context(system, platform):
    """
    Create an OpenMM Context for repeated single-point evaluations.
    """
    integrator = LangevinMiddleIntegrator(
        TEMPERATURE,
        FRICTION,
        TIMESTEP,
    )

    context = Context(system, integrator, platform)
    return context, integrator


def context_energy_kcal(context, positions_nm):
    """
    Evaluate potential energy in kcal/mol.
    """
    context.setPositions(
        positions_nm * unit.nanometer
    )

    state = context.getState(getEnergy=True)

    return state.getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole
    )


def make_energy_contexts(structure):
    """
    Create contexts for:
      - total MM energy
      - bonded-only
      - electrostatic + bonded
      - vdW + bonded
      - GBSA total
    """
    platform = choose_platform()

    vacuum = create_vacuum_system(structure)

    bonded_system = modify_nonbonded(vacuum, "bonded")
    electrostatic_system = modify_nonbonded(vacuum, "electrostatic")
    vdw_system = modify_nonbonded(vacuum, "vdw")

    gbsa_system = create_gbsa_system(structure)

    systems = {
        "mm": vacuum,
        "bonded": bonded_system,
        "electrostatic": electrostatic_system,
        "vdw": vdw_system,
        "gbsa": gbsa_system,
    }

    contexts = {}
    integrators = []

    for key, system in systems.items():
        ctx, integrator = make_context(system, platform)
        contexts[key] = ctx
        integrators.append(integrator)

    return contexts, integrators


def evaluate_species_energy(structure, positions_nm):
    """
    Evaluate all MM/GBSA components for one species at one snapshot.

    Returns:
        E_MM
        E_internal
        E_electrostatic
        E_vdw
        G_solv
        G_total
    """
    contexts, integrators = make_energy_contexts(structure)

    try:
        e_mm = context_energy_kcal(
            contexts["mm"],
            positions_nm,
        )

        e_internal = context_energy_kcal(
            contexts["bonded"],
            positions_nm,
        )

        e_electrostatic_plus_internal = context_energy_kcal(
            contexts["electrostatic"],
            positions_nm,
        )

        e_vdw_plus_internal = context_energy_kcal(
            contexts["vdw"],
            positions_nm,
        )

        e_gbsa_total = context_energy_kcal(
            contexts["gbsa"],
            positions_nm,
        )

    finally:
        for integrator in integrators:
            del integrator

    e_electrostatic = (
        e_electrostatic_plus_internal - e_internal
    )

    e_vdw = (
        e_vdw_plus_internal - e_internal
    )

    g_solv = e_gbsa_total - e_mm

    return {
        "E_MM": e_mm,
        "E_internal": e_internal,
        "E_electrostatic": e_electrostatic,
        "E_vdw": e_vdw,
        "G_solv": g_solv,
        "G_total": e_gbsa_total,
    }


def strip_solute_structures(complex_parm):
    """
    Generate receptor, ligand and complex structures from the parametrized
    complex.

    Water and ions are removed before MM/GBSA.
    """
    ligand_resname = LIGAND_ID[:3].upper()

    ligand_mask = np.array([
        atom.residue.name.strip().upper() == ligand_resname
        for atom in complex_parm.atoms
    ])

    if not ligand_mask.any():
        raise RuntimeError("Ligand could not be identified in ParmEd structure.")

    # Protein atoms only.
    receptor_mask = np.array([
        atom.residue.name.upper() not in {
            ligand_resname,
            "WAT",
            "HOH",
            "SOL",
            "TIP3",
            "NA",
            "CL",
        }
        and atom.residue.name.upper() not in {"NA+", "CL-"}
        for atom in complex_parm.atoms
    ])

    # More robust receptor selection: retain atoms that are not water/ions and are not ligand. This intentionally includes protein heteroatoms that are part of the parametrized receptor.
    solvent_names = {
        "WAT",
        "HOH",
        "SOL",
        "TIP3",
        "TIP3P",
        "NA",
        "CL",
        "SOD",
        "CLA",
    }

    receptor_mask = np.array([
        atom.residue.name.strip().upper() not in solvent_names
        and atom.residue.name.strip().upper() != ligand_resname
        for atom in complex_parm.atoms
    ])

    ligand = complex_parm[ligand_mask]
    receptor = complex_parm[receptor_mask]

    # MM/GBSA complex contains receptor + ligand but no bulk solvent/ions.
    complex_solute_mask = receptor_mask | ligand_mask
    complex_solute = complex_parm[complex_solute_mask]

    return complex_solute, receptor, ligand


def run_mmgbsa(analysis_traj, complex_parm):
    """
    Perform single-trajectory MM/GBSA over the selected MDTraj frames.
    """
    complex_solute, receptor, ligand = strip_solute_structures(
        complex_parm
    )

    (
        protein_atoms,
        ligand_atoms,
        ligand_heavy,
        protein_backbone,
        protein_ca,
    ) = identify_indices(analysis_traj)

    ligand_resname = LIGAND_ID[:3].upper()

    # Build MDTraj atom maps from the full trajectory to the stripped
    # complex/receptor/ligand structures.
    full_solute_indices = np.concatenate(
        [protein_atoms, ligand_atoms]
    )

    # The ParmEd complex order is receptor first and ligand last. The MD
    # trajectory has the same order because it was generated from complex_parm.
    #
    # We therefore construct a map directly from the ParmEd structure order.
    ligand_parmed_indices = [
        i for i, atom in enumerate(complex_parm.atoms)
        if atom.residue.name.strip().upper() == ligand_resname
    ]

    receptor_parmed_indices = [
        i for i, atom in enumerate(complex_parm.atoms)
        if atom.residue.name.strip().upper() not in {
            ligand_resname,
            "WAT",
            "HOH",
            "SOL",
            "TIP3",
            "TIP3P",
            "NA",
            "CL",
            "SOD",
            "CLA",
        }
    ]

    solute_parmed_indices = (
        receptor_parmed_indices + ligand_parmed_indices
    )

    # Convert ParmEd/Å coordinates into nm and use exactly the same order as
    # the trajectory.
    frame_indices = np.arange(analysis_traj.n_frames)

    rows = []

    print(
        f"Running MM/GBSA on {len(frame_indices)} frames "
        f"(stride={MMGBSA_FRAME_STRIDE})..."
    )

    for n, frame in enumerate(frame_indices, start=1):
        xyz_full_nm = analysis_traj.xyz[frame]

        complex_positions = xyz_full_nm[solute_parmed_indices]
        receptor_positions = xyz_full_nm[receptor_parmed_indices]
        ligand_positions = xyz_full_nm[ligand_parmed_indices]

        c = evaluate_species_energy(
            complex_solute,
            complex_positions,
        )

        r = evaluate_species_energy(
            receptor,
            receptor_positions,
        )

        l = evaluate_species_energy(
            ligand,
            ligand_positions,
        )

        row = {
            "frame": int(frame),
            "time_ns": float(analysis_traj.time[frame] / 1000.0),

            "E_MM_complex": c["E_MM"],
            "E_MM_receptor": r["E_MM"],
            "E_MM_ligand": l["E_MM"],

            "E_internal_complex": c["E_internal"],
            "E_internal_receptor": r["E_internal"],
            "E_internal_ligand": l["E_internal"],

            "E_electrostatic_complex": c["E_electrostatic"],
            "E_electrostatic_receptor": r["E_electrostatic"],
            "E_electrostatic_ligand": l["E_electrostatic"],

            "E_vdw_complex": c["E_vdw"],
            "E_vdw_receptor": r["E_vdw"],
            "E_vdw_ligand": l["E_vdw"],

            "G_solv_complex": c["G_solv"],
            "G_solv_receptor": r["G_solv"],
            "G_solv_ligand": l["G_solv"],

            "G_total_complex": c["G_total"],
            "G_total_receptor": r["G_total"],
            "G_total_ligand": l["G_total"],
        }


        # Single-trajectory MM/GBSA thermodynamic cycle:
        #
        # ΔGbind = Gcomplex - Greceptor - Gligand
        #
        # No -TΔS term is included.


        row["dE_MM"] = (
            row["E_MM_complex"]
            - row["E_MM_receptor"]
            - row["E_MM_ligand"]
        )

        row["dE_internal"] = (
            row["E_internal_complex"]
            - row["E_internal_receptor"]
            - row["E_internal_ligand"]
        )

        row["dE_electrostatic"] = (
            row["E_electrostatic_complex"]
            - row["E_electrostatic_receptor"]
            - row["E_electrostatic_ligand"]
        )

        row["dE_vdw"] = (
            row["E_vdw_complex"]
            - row["E_vdw_receptor"]
            - row["E_vdw_ligand"]
        )

        row["dG_solv"] = (
            row["G_solv_complex"]
            - row["G_solv_receptor"]
            - row["G_solv_ligand"]
        )

        row["dG_MMGBSA"] = (
            row["G_total_complex"]
            - row["G_total_receptor"]
            - row["G_total_ligand"]
        )

        rows.append(row)

        if n % 10 == 0 or n == len(frame_indices):
            print(f"  {n}/{len(frame_indices)}")

    df = pd.DataFrame(rows)
    df.to_csv(ENERGY_CSV, index=False)

    return df


def summarize_mmgbsa(df):
    """
    Calculate mean, SD, SEM and cumulative mean for the MM/GBSA ensemble.
    """
    components = [
        "dE_MM",
        "dE_internal",
        "dE_electrostatic",
        "dE_vdw",
        "dG_solv",
        "dG_MMGBSA",
    ]

    rows = []

    for component in components:
        values = df[component].to_numpy(dtype=float)

        rows.append({
            "component": component,
            "mean_kcal_mol": np.mean(values),
            "std_kcal_mol": np.std(values, ddof=1),
            "sem_kcal_mol": (
                np.std(values, ddof=1) / np.sqrt(len(values))
            ),
            "n_frames": len(values),
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(SUMMARY_CSV, index=False)

    return summary


def plot_mmgbsa(df):
    """
    Plot per-frame MM/GBSA components and convergence.
    """
    time = df["time_ns"].to_numpy()


    # Per-frame energetic components.

    plt.figure(figsize=(10, 6))
    plt.plot(time, df["dE_vdw"], label="ΔE_vdW")
    plt.plot(time, df["dE_electrostatic"], label="ΔE_electrostatic")
    plt.plot(time, df["dG_solv"], label="ΔG_solv")
    plt.plot(time, df["dG_MMGBSA"], label="ΔG_MMGBSA")
    plt.xlabel("Time (ns)")
    plt.ylabel("Energy (kcal/mol)")
    plt.title("MM/GBSA energy components")
    plt.legend()
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mmgbsa_components.png", dpi=300)
    plt.close()


    # Cumulative mean of ΔG_MMGBSA.

    cumulative = df["dG_MMGBSA"].expanding().mean()

    plt.figure(figsize=(8, 5))
    plt.plot(time, cumulative)
    plt.axhline(
        cumulative.iloc[-1],
        linestyle="--",
        linewidth=1,
    )
    plt.xlabel("Time (ns)")
    plt.ylabel("Cumulative mean ΔG_MMGBSA (kcal/mol)")
    plt.title("MM/GBSA convergence")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mmgbsa_convergence.png", dpi=300)
    plt.close()


    # Distribution of ΔG_MMGBSA.

    plt.figure(figsize=(8, 5))
    plt.hist(
        df["dG_MMGBSA"],
        bins=30,
        density=True,
    )
    plt.xlabel("ΔG_MMGBSA (kcal/mol)")
    plt.ylabel("Density")
    plt.title("MM/GBSA energy distribution")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mmgbsa_distribution.png", dpi=300)
    plt.close()



# MAIN WORKFLOW


print("=" * 72)
print("START MM/GBSA WORKFLOW")
print("=" * 72)

print("\n[1/8] Ligand preparation")
ligand_top = create_ligand_topology()

print("\n[2/8] Selecting crystallographic pocket waters")
pocket_water = prepare_pocket_waters()

print("\n[3/8] Building dry complex")
build_dry_complex(pocket_water)

print("\n[4/8] Adding hydrogens, water and ions")
modeller, forcefield = solvate_with_openmm(COMPLEX_DRY_PDB)

print("\n[5/8] Building parametrized complex")
complex_parm = parameterize_explicit_system(
    modeller,
    forcefield,
    ligand_top,
)

print(
    f"Parametrized complex: {len(complex_parm.atoms)} atoms, "
    f"{len(complex_parm.residues)} residues"
)

print("\n[6/8] Running OpenMM MD")
run_md(complex_parm)

print("\n[7/8] MDTraj analysis")
trajectory_topology = build_analysis_topology(complex_parm)

full_traj, analysis_traj = load_production_trajectory(
    trajectory_topology
)

indices = identify_indices(analysis_traj)

trajectory_metrics = plot_trajectory_analysis(
    full_traj,
    analysis_traj,
    indices,
)

print("\n[8/8] MM/GBSA")
mmgbsa_df = run_mmgbsa(
    analysis_traj,
    complex_parm,
)

summary = summarize_mmgbsa(mmgbsa_df)
plot_mmgbsa(mmgbsa_df)

print("\n" + "=" * 72)
print("MM/GBSA SUMMARY")
print("=" * 72)
print(summary.to_string(index=False))

print("\nFiles written to:")
print(RESULTS_DIR)
print("\nDONE.")
