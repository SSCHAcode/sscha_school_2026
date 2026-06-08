import warnings

import ase
import ase.io
import numpy as np
from ase import Atoms
from ase.eos import calculate_eos
from ase.filters import ExpCellFilter
from ase.optimize import QuasiNewton
from ase.units import kJ
from quippy.potential import Potential

# Import the CellConstructor and SSCHA modules used throughout this tutorial.
import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.Structure
import sscha
import sscha.Ensemble
import sscha.Relax
import sscha.SchaMinimizer

# NumPy moved ComplexWarning in newer versions. This block keeps the script
# compatible with both older and newer NumPy releases.
try:
    ComplexWarning = np.exceptions.ComplexWarning
except AttributeError:
    ComplexWarning = np.ComplexWarning

warnings.filterwarnings("ignore", category=ComplexWarning)

def get_starting_dynamical_matrices(structure_filename, potential, supercell):
    """Relax the input structure and compute an initial dynamical matrix.

    Parameters
    ----------
    structure_filename : str
        Path to the input structure in Quantum ESPRESSO format.
    potential : quippy.potential.Potential
        Interatomic potential used to compute energies and forces.
    supercell : list[int]
        Supercell used for the finite-displacement phonon calculation.

    Returns
    -------
    dyn : cellconstructor.Phonons.Phonons
        Initial auxiliary dynamical matrix.
    bulk : float
        Static bulk modulus, used during variable-cell SSCHA relaxation.
    """
    # Read the input structure with CellConstructor.
    structure = CC.Structure.Structure()
    structure.read_scf(structure_filename)

    # Convert the CellConstructor structure to an ASE Atoms object so that
    # ASE optimizers and the QUIP/GAP potential can be used.
    atoms = Atoms(
        cell=structure.unit_cell,
        positions=structure.coords,
        symbols=structure.atoms,
        pbc=True,
    )
    atoms.set_calculator(potential)

    # Relax both atomic positions and cell vectors.
    ecf = ExpCellFilter(atoms)
    qn = QuasiNewton(ecf)
    qn.run(fmax=0.0005)

    # Convert the relaxed ASE structure back to CellConstructor format.
    structure = CC.Structure.Structure()
    structure.generate_from_ase_atoms(atoms, get_masses=True)

    # Compute an initial dynamical matrix by finite displacements.
    dyn = CC.Phonons.compute_phonons_finite_displacements(
        structure,
        potential,
        supercell=supercell,
    )
    dyn.Symmetrize()

    # Remove imaginary frequencies from the starting guess. This produces a
    # positive-definite auxiliary dynamical matrix for the SSCHA minimization.
    dyn.ForcePositiveDefinite()

    # Estimate the static bulk modulus from an equation-of-state fit. This is
    # needed by the variable-cell SSCHA relaxation below.
    bulk = get_bulk_modulus(dyn, potential)

    return dyn, bulk

def get_bulk_modulus(dyn, pot):
    atoms = dyn.structure.get_ase_atoms()
    atoms.set_calculator(pot)
    eos = calculate_eos(atoms)
    v0, e0, B = eos.fit()
    return B / kJ * 1.0e24

# Input files and calculation settings. Adjust these paths for your machine.
src = "/home/dangic/Documents/Posters_and_talks/SSCHA_school_2026/sscha_school_2026/Materials/"
structure_name = "CsPbI3_cubic_phase.scf"
pot = Potential(param_filename=src + "GAP_1.xml", calc_args="local_gap_variance")
supercell = [2, 2, 2]
dyn_filename = 'sscha_auxiliary_dyn_'
# SSCHA relaxation settings.
temperature = 300.0
nconf = 500
max_pop = 1000
relax = False
numconf = 2000

# Generate the initial auxiliary dynamical matrix and the bulk modulus.
if(dyn_filename is None):
    dyn, bulk = get_starting_dynamical_matrices(src + structure_name, pot, supercell)
else:
    dyn = CC.Phonons.Phonons(src + dyn_filename, nqirr = 4)
    bulk = get_bulk_modulus(dyn, pot)

if(relax):
    # Generate the initial SSCHA ensemble from the starting dynamical matrix.
    ensemble = sscha.Ensemble.Ensemble(
        dyn,
        T0=temperature,
        supercell=dyn.GetSupercell(),
    )
    ensemble.generate(N=nconf)

    # Set up the SSCHA minimizer.
    minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    minimizer.min_step_dyn = 0.001
    minimizer.kong_liu_ratio = 0.5
    minimizer.meaningful_factor = 0.001
    minimizer.max_ka = 100000

    # Perform the variable-cell SSCHA relaxation.
    relax = sscha.Relax.SSCHA(
        minimizer,
        ase_calculator=pot,
        N_configs=nconf,
        max_pop=max_pop,
        save_ensemble=True,
    )
    relax.vc_relax(
        static_bulk_modulus=bulk,
        ensemble_loc="directory_of_the_ensemble",
    )

    dyn = relax.minim.dyn.Copy()


# Generate a new ensemble using the relaxed SSCHA dynamical matrix.
new_ensemble = sscha.Ensemble.Ensemble(
    dyn,
    T0=temperature,
    supercell=dyn.GetSupercell(),
)
new_ensemble.generate(N=numconf)

# Compute energies, forces, and stresses for the generated configurations.
new_ensemble.compute_ensemble(
    pot,
    compute_stress=True,
    stress_numerical=False,
    cluster=None,
    verbose=True,
)

# Minimize the free energy for the new ensemble. We keep the structure fixed
# here because the variable-cell relaxation was already performed above.
new_minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(new_ensemble)
new_minimizer.minim_struct = False
new_minimizer.set_minimization_step(0.001)
new_minimizer.meaningful_factor = 0.001
new_minimizer.max_ka = 10000
new_minimizer.init()
new_minimizer.run()

# Save the converged auxiliary dynamical matrix in Quantum ESPRESSO format.
new_minimizer.dyn.save_qe("final_dyn_" + str(numconf) + "_")

# Reweight the ensemble with the updated dynamical matrix.
new_ensemble.update_weights(new_minimizer.dyn, temperature)

# Calculate the free-energy Hessian and the third-order tensor.
# Setting return_d3=True returns the third-order force constants needed for
# the thermal-conductivity calculation.
dyn_hessian, d3_tensor = new_ensemble.get_free_energy_hessian(
    include_v4=False,
    get_full_hessian=True,
    return_d3=True,
)

# Save the third-order tensor and the Hessian for later use.
np.save("d3_" + str(numconf) + ".npy", d3_tensor)
dyn_hessian.save_qe("hessian_dyn_" + str(numconf) + "_")
