# Hands-on Session 9: Thermal Conductivity Calculations with the SSCHA

In previous lessons, we saw how to calculate vibrational properties of materials using the stochastic self-consistent harmonic approximation (SSCHA). In this tutorial, we will use that knowledge to calculate the lattice thermal conductivity of materials.

To do this, we need the SSCHA dynamical matrices (**auxiliary dynamical matrices, not free-energy Hessians**) and the third-order force constants. These quantities allow us to calculate the harmonic properties, such as phonon frequencies and group velocities, and the anharmonic properties, such as phonon lifetimes and spectral functions. Together, these are the ingredients required for lattice thermal-conductivity calculations.

In this tutorial, we will:

- Use SSCHA dynamical matrices and third-order force constants to calculate lattice thermal conductivity with the Green-Kubo (GK) and single relaxation time approximation (SRTA) approaches.
- Learn how to converge lattice thermal-conductivity calculations.
- Manipulate `ThermalConductivity` objects to analyze the calculated transport properties.

The basic workflow is:

1. Relax the crystal structure and obtain SSCHA dynamical matrices in the temperature range of interest.
2. Use the relaxed SSCHA dynamical matrices to calculate the SSCHA free-energy Hessian, confirm dynamical stability, and obtain third-order force constants.
3. Calculate the lattice thermal conductivity and check its convergence with respect to:
   - the **q**-point grid used for the thermal-conductivity calculation;
   - the size of the SSCHA supercell;
   - the number of configurations used to calculate the third-order force constants;
   - the smearing parameter used in the phonon self-energy calculation, which is closely related to the size of the **q**-point grid when using an *adaptive* smearing scheme;
   - the density of frequency sampling controlled by the `ne` parameter when using the Green-Kubo method.

## Lattice Thermal Conductivity of CsPbI<sub>3</sub>

We first calculate the second- and third-order force constants.

```python
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
    eos = calculate_eos(atoms)
    v0, e0, B = eos.fit()
    bulk = B / kJ * 1.0e24

    return dyn, bulk


# Input files and calculation settings. Adjust these paths for your machine.
src = "/scratch/ddangic/sscha_school_2026/Materials/"
structure_name = "CsPbI3_cubic_phase.scf"
pot = Potential(param_filename=src + "GAP_1.xml", calc_args="local_gap_variance")
supercell = [2, 2, 2]

# Generate the initial auxiliary dynamical matrix and the bulk modulus.
dyn, bulk = get_starting_dynamical_matrices(src + structure_name, pot, supercell)

# SSCHA relaxation settings.
temperature = 300.0
nconf = 500
max_pop = 1000

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

# Recompute the SSCHA ensemble with an increasing number of configurations.
# This loop is used to test the convergence of the third-order force constants
# and the corresponding thermal-conductivity results.
for i in range(10):
    numconf = (i + 1) * 1000

    # Generate a new ensemble using the relaxed SSCHA dynamical matrix.
    new_ensemble = sscha.Ensemble.Ensemble(
        relax.minim.dyn,
        T0=temperature,
        supercell=relax.minim.dyn.GetSupercell(),
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
```

In this example, we used a 2 × 2 × 2 supercell to calculate the second- and third-order force constants. The third-order force constants were calculated using an increasing number of configurations in order to test the convergence of the ensemble average.

The second- and third-order force constants must be converged with respect to both the supercell size and the number of configurations. Figure 1 shows the results for harmonic properties. The phonon band structure interpolated from the 2 × 2 × 2 supercell is converged with respect to the number of configurations at around 2000 configurations. However, the results also need to be converged with respect to the supercell size, and they do not appear to be converged even for the 4 × 4 × 4 supercell.

<figure>
  <img src="./phonons.jpg" alt="Phonon band structures calculated with different supercell sizes and numbers of SSCHA configurations.">
  <figcaption>Figure 1. Convergence study of the phonon band structure with respect to supercell size and the number of configurations.</figcaption>
</figure>

First, we calculate the lattice thermal conductivity using the single relaxation time approximation (**SRTA**). In this approximation, the thermal-conductivity tensor is written as

```math
\kappa^{xy}
=
\frac{1}{N V}
\sum_{\mathbf{q}, j}
v^{x}_{\mathbf{q}, j}
v^{y}_{\mathbf{q}, j}
c_{\mathbf{q}, j}
\tau_{\mathbf{q}, j}.
```

Here, $`N`$ is the number of **q** points,  $`V`$ is the unit-cell volume,  $`v^{x}_{\mathbf{q}, j}`$ and  $`v^{y}_{\mathbf{q}, j}`$ are the phonon group-velocity components,  $`c_{\mathbf{q}, j}`$ is the mode heat capacity, and  $`\tau_{\mathbf{q}, j}`$ is the phonon lifetime for phonon mode  $`j`$ at wave vector  $`\mathbf{q}`$. To calculate the thermal conductivity using force constants obtained with different numbers of configurations, we can use the following script:

```python
from __future__ import division, print_function

import time

import numpy as np

import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Phonons
import cellconstructor.ThermalConductivity

# Number of irreducible q points/files in the saved dynamical matrix.
nqirr = 4

# Useful conversion factors from CellConstructor.
SSCHA_TO_MS = cellconstructor.ThermalConductivity.SSCHA_TO_MS
RY_TO_THZ = cellconstructor.ThermalConductivity.SSCHA_TO_THZ

# Use the same grid for phonon properties and scattering processes in this
# example. In production calculations, these two grids should be converged
# separately.
meshnum = 10
mesh = [meshnum, meshnum, meshnum]

# Loop over force constants obtained from different ensemble sizes. This starts
# from 2000 configurations because range(1, 10) gives i = 1, ..., 9.
for i in range(1, 10):
    numconf = (i + 1) * 1000

    # Read the auxiliary SSCHA dynamical matrix saved in the previous step.
    dyn_prefix = "final_dyn_" + str(numconf) + "_"
    dyn = CC.Phonons.Phonons(dyn_prefix, nqirr)
    supercell = dyn.GetSupercell()

    # Create the third-order force-constant object using the same supercell
    # used to generate the SSCHA dynamical matrix.
    fc3 = CC.ForceTensor.Tensor3(
        dyn.structure,
        dyn.structure.generate_supercell(supercell),
        supercell,
    )

    # Load the third-order tensor and use it to initialize fc3.
    d3_name = "d3_" + str(numconf) + ".npy"
    d3 = np.load(d3_name)
    fc3.SetupFromTensor(d3)

    # Center the third-order force constants. This step is necessary 
    # in order to interpolate 3rd order force constants.
    fc3 = CC.ThermalConductivity.centering_fc3(fc3)

    # Define the ThermalConductivity object. The kpoint_grid is the grid on which
    # phonon quantities are calculated on, while scattering_grid controls 
    # the sampling of scattering processes.
    tc = CC.ThermalConductivity.ThermalConductivity(
        dyn,
        fc3,
        kpoint_grid=mesh,
        scattering_grid=mesh,
        smearing_scale=None,
        smearing_type="adaptive",
        cp_mode="quantum",
        off_diag=True,
    )

    # Calculate the thermal conductivity at 300 K.
    temperatures = np.linspace(300, 300, 1, dtype=float)
    start_time = time.time()

    # Compute harmonic quantities such as phonon frequencies and group
    # velocities, then write them to text files.
    tc.setup_harmonic_properties()
    tc.write_harmonic_properties_to_file()

    # Calculate the lattice thermal conductivity in the SRTA approximation.
    # Lifetimes are written to file because write_lifetimes=True.
    tc.calculate_kappa(
        mode="SRTA",
        temperatures=temperatures,
        write_lifetimes=True,
        gauss_smearing=True,
        offdiag_mode="perturbative",
        kappa_filename="ConfThermal_conductivity_SRTA_" + str(numconf),
        lf_method="fortran-LA",
    )
```

The most important parts of the script are:

- The mesh used to calculate phonon properties is defined to be the same as the mesh used to calculate scattering processes. This is controlled by the variable `mesh`. These grids do not have to be identical. In many cases, `scattering_grid` can be much smaller than `kpoint_grid`. **Both grids must be converged.**
- We use a smearing approach to satisfy energy-conservation laws. There are two methods: constant and adaptive. For `smearing_type = "constant"`, the smearing value must be provided in **Ry** as an argument to `setup_harmonic_properties`. For adaptive smearing, as used here, the smearing value is different for different phonon modes and is based on the phonon density of states and **q**-point grid density. A global scaling factor can still be provided through `smearing_scale`; `smearing_scale = 1.0` often works well. **The smearing parameters must also be converged.**
- The `off_diag` variable controls whether the calculation includes what is often called *coherent transport*. This contribution can be important for strongly anharmonic materials with significant bunching of phonon modes.
- The function `calculate_kappa` performs most of the transport calculation. Its main options are:
  - `mode` defines the method used to calculate the lattice thermal conductivity. The options are **SRTA**, which corresponds to the Boltzmann transport equation in the single relaxation time approximation, and **GK** ([Dangic et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.104314)), which is the Green-Kubo method based on phonon spectral functions rather than phonon lifetimes. These two modes should give similar results in weakly anharmonic materials, but they can differ in strongly anharmonic materials.
  - `gauss_smearing` defines how energy conservation is treated in the self-energy calculation. If `True`, Gaussian functions are used. If `False`, Lorentzian functions are used. With Gaussian smearing, the real part of the self-energy is calculated using the Kramers-Kronig transformation.
  - `offdiag_mode` defines how coherent transport is calculated when `mode = "SRTA"`. The available options are **wigner** ([Simoncelli et al.](https://www.nature.com/articles/s41567-019-0520-x)), **gk** ([Isaeva et al.](https://www.nature.com/articles/s41467-019-11572-4)), and **perturbative** ([Dangic et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.104314)). If `mode = "GK"`, coherent transport is included naturally.
  - `lf_method` defines how phonon lifetimes are calculated when `mode = "SRTA"`. In most cases, keep the `fortran-` prefix and add either `LA` or `P`; these options should give similar results. Another option is `SC`, where phonon lifetimes are solved self-consistently and phonon lineshifts are included.
  - `ne` defines the number of frequency steps used when calculating phonon lineshapes. It is also important for `lf_method = "SC"`, because the self-consistent equation is solved on a frequency grid and the real and imaginary parts of the self-energy are linearly interpolated. Larger values are generally more accurate. **The results must be converged with respect to `ne`.**

In this example, we use the *SRTA* approximation, which is the weak-anharmonicity limit of the Green-Kubo equation. We use SSCHA auxiliary dynamical matrices because they form the basis for the dynamical extension of the SSCHA, derived from time-dependent SSCHA, which gives access to phonon lifetimes and lineshapes.

By setting `off_diag=True`, we include the coherent contribution to the lattice thermal conductivity. This contribution is important in strongly anharmonic materials, alloys, and amorphous solids. The calculation should take a few minutes, and the results are saved using the filename prefix specified by `kappa_filename`.

The thermal conductivity should also be calculated using different **q**-point grids. Figure 2 includes results obtained for a 3 × 3 × 3 grid. These calculations do not appear to be converged at the 3 × 3 × 3 level, so a larger supercell would be required for a production calculation.

However, third-order force constants often converge faster with respect to supercell size than second-order force constants. Therefore, one can sometimes combine dynamical matrices computed on a larger supercell with third-order force constants computed on a smaller supercell. This approach is illustrated in the script below:
```python
from __future__ import print_function
from __future__ import division

import time

import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.ThermalConductivity

SSCHA_TO_MS = cellconstructor.ThermalConductivity.SSCHA_TO_MS
RY_TO_THZ = cellconstructor.ThermalConductivity.SSCHA_TO_THZ

# Define the q-point grid used for both harmonic properties and scattering processes.
meshnum = 10
mesh = [meshnum, meshnum, meshnum]

# Define the temperature at which the thermal conductivity will be calculated.
temperatures = np.linspace(300, 300, 1, dtype=float)

# Read the dynamical matrix corresponding to the same supercell used to compute
# the third-order force constants.
dyn_prefix = "final_dyn_2x2x2_"
nqirr = 4
dyn_for_fc3 = CC.Phonons.Phonons(dyn_prefix, nqirr)

# Initialize and load the third-order force-constant tensor.
supercell = dyn_for_fc3.GetSupercell()
fc3 = CC.ForceTensor.Tensor3(
    dyn_for_fc3.structure,
    dyn_for_fc3.structure.generate_supercell(supercell),
    supercell,
)

d3_name = "d3_2x2x2.npy"
d3 = np.load(d3_name)
fc3.SetupFromTensor(d3)

# Center the third-order force constants.
fc3 = CC.ThermalConductivity.centering_fc3(fc3)

# Read a dynamical matrix computed using a larger supercell.
# This matrix will provide the harmonic properties used in the transport calculation.
dyn_prefix = "./final_dyn_4x4x4_"
nqirr = 10
dyn = CC.Phonons.Phonons(dyn_prefix, nqirr)

# Build the ThermalConductivity object using the large-supercell dynamical matrix
# and the small-supercell third-order force constants.
tc = CC.ThermalConductivity.ThermalConductivity(
    dyn,
    fc3,
    kpoint_grid=mesh,
    scattering_grid=mesh,
    smearing_scale=None,
    smearing_type="adaptive",
    cp_mode="quantum",
    off_diag=True,
)

start_time = time.time()

# Compute harmonic properties and write them to file.
tc.setup_harmonic_properties()
tc.write_harmonic_properties_to_file()

# Calculate the thermal conductivity in the SRTA approximation.
tc.calculate_kappa(
    mode="SRTA",
    temperatures=temperatures,
    write_lifetimes=True,
    gauss_smearing=True,
    offdiag_mode="perturbative",
    kappa_filename="Thermal_conductivity_SRTA",
    lf_method="fortran-LA",
)
```

<figure>
  <img src="./tcs.jpg" alt="Thermal conductivity convergence with respect to q-point grid, supercell size, and number of configurations.">
  <figcaption>Figure 2. Convergence study of thermal conductivity with respect to q-point sampling, supercell size, and number of configurations.</figcaption>
</figure>

Finally, we can calculate the lattice thermal conductivity with the Green-Kubo method:

```python
from __future__ import division, print_function

import time

import numpy as np

import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Phonons
import cellconstructor.ThermalConductivity

# Use the dynamical matrix and third-order force constants obtained with
# 10000 configurations.
dyn_prefix = "final_dyn_10000_"
nqirr = 4

# Useful conversion factors.
SSCHA_TO_MS = cellconstructor.ThermalConductivity.SSCHA_TO_MS
RY_TO_THZ = cellconstructor.ThermalConductivity.SSCHA_TO_THZ

# Read the auxiliary SSCHA dynamical matrix.
dyn = CC.Phonons.Phonons(dyn_prefix, nqirr)
supercell = dyn.GetSupercell()

# Build and initialize the third-order force-constant tensor.
fc3 = CC.ForceTensor.Tensor3(
    dyn.structure,
    dyn.structure.generate_supercell(supercell),
    supercell,
)
d3 = np.load("d3_10000.npy")
fc3.SetupFromTensor(d3)
fc3 = CC.ThermalConductivity.centering_fc3(fc3)

# Define the q-point grid for the Green-Kubo calculation.
meshnum = 10
mesh = [meshnum, meshnum, meshnum]

# Set up the ThermalConductivity object. With mode="GK", the method uses
# phonon spectral functions rather than phonon lifetimes.
tc = CC.ThermalConductivity.ThermalConductivity(
    dyn,
    fc3,
    kpoint_grid=mesh,
    scattering_grid=mesh,
    smearing_scale=1.0,
    smearing_type="adaptive",
    cp_mode="quantum",
    off_diag=True,
)

# Compute harmonic quantities once. They are reused for all values of ne below.
temperatures = np.linspace(300, 300, 1, dtype=float)
start_time = time.time()
tc.setup_harmonic_properties()
tc.write_harmonic_properties_to_file()

# Converge the Green-Kubo calculation with respect to the number of frequency
# points used to integrate overlaps of phonon spectral functions.
for i in range(5):
    ne = i * 2000 + 1000

    tc.calculate_kappa(
        mode="GK",
        temperatures=temperatures,
        write_lifetimes=False,
        gauss_smearing=False,
        ne=ne,
        offdiag_mode="perturbative",
        kappa_filename="Thermal_conductivity_GK_" + str(ne),
        lf_method="fortran-LA",
    )

# Save the ThermalConductivity object for post-processing.
tc.save_pickle()
```

The Green-Kubo calculation introduces an additional convergence parameter: the number of frequency points used to integrate the overlap of phonon spectral functions. The frequency spacing is approximately $2.1 \omega_D / n_e$, where $\omega_D$ is the highest frequency on the **q**-point grid. For an accurate estimate of the integrals, this frequency spacing should be smaller than the smallest phonon linewidth. Therefore, materials with higher thermal conductivity will likely require larger values of `ne`.

The SRTA and GK methods usually agree well for weakly anharmonic materials. For strongly anharmonic materials, however, their results can differ significantly.

The final line of the previous script saves the `ThermalConductivity` object, which can later be used to analyze the calculated transport quantities. The snippet below reloads this object and writes the phonon spectral function at $\Gamma$:

```python
from __future__ import division, print_function

import time

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Phonons
import cellconstructor.ThermalConductivity

# Reload the ThermalConductivity object saved with tc.save_pickle().
tc = CC.ThermalConductivity.load_thermal_conductivity()

# Temperature-dependent quantities, such as lineshapes, are stored in
# dictionaries. The keys are strings formatted as format(T, ".1f").
keys = list(tc.lineshapes.keys())

# Find the Gamma point in the q-point grid.
for iqpt in range(tc.nkpt):
    if np.linalg.norm(tc.k_points[iqpt]) == 0.0:
        break

# Build the frequency grid associated with the saved lineshape data.
energies = (
    np.arange(len(tc.lineshapes[keys[-1]][0, 0]), dtype=float) * tc.delta_omega
    + tc.delta_omega
)

# Write the phonon spectral function at Gamma to a text file.
tc.write_lineshape(
    "Lineshape_at_Gamma",
    tc.lineshapes[keys[-1]][iqpt],
    iqpt,
    energies,
    "no",
)
```

Temperature-dependent transport properties, such as phonon lifetimes, heat capacities, and phonon lineshapes, are stored as dictionaries with keys corresponding to the temperature `T` formatted as `format(T, ".1f")`. The example above writes the phonon spectral function at $\Gamma$.

Other relevant quantities, such as phonon frequencies and group velocities, can also be accessed from the `ThermalConductivity` object. In addition, the object provides built-in functions for calculating the phonon density of states, phonon spectral conductivity, and phonon Grüneisen parameters in cubic systems

## Exercise

Estimate the lattice thermal conductivity of CsPbI<sub>3</sub> in the cubic $Pm\bar{3}m$ phase.

Before computing the thermal conductivity, first check the dynamical stability of this phase. In particular, determine the temperature range in which the cubic phase is dynamically stable and therefore the thermal-conductivity calculation is physically meaningful.

Then compare the lattice thermal conductivity obtained using:

- the single relaxation time approximation (**SRTA**), and
- the Green-Kubo (**GK**) approach.

Discuss the following points:

1. At what temperature does the cubic $Pm\bar{3}m$ phase become dynamically stable?
2. How does the thermal conductivity depend on temperature?
3. Do the SRTA and GK results agree? If not, which method gives a larger thermal conductivity and why?
4. What does the comparison between SRTA and GK tell you about the degree of anharmonicity in cubic CsPbI<sub>3</sub>?
