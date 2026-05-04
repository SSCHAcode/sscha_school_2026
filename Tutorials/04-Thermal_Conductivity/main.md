# Hands-on Session 9 - Thermal Conductivity Calculations with the SSCHA

In previous lessons, we saw how to calculate vibrational properties of materials using the SSCHA. Now we will use this acquired knowledge to calculate the lattice thermal conductivity of materials. We will need dynamical matrices (**auxiliary ones, not Hessians**) and the third-order force constants (we already calculated them when we checked the dynamical stability of the system). With these, we can calculate the harmonic (phonon frequencies and phonon group velocities) and anharmonic properties (phonon lifetimes and spectral functions) of materials, which is all we need to calculate lattice thermal conductivity.

## Lattice Thermal Conductivity of Silicon

As a first exercise, let us calculate the lattice thermal conductivity of silicon. Silicon is a very harmonic material, which means its lattice thermal conductivity is very high. This also makes it a good test case to check the equivalence of the Green-Kubo and Boltzmann transport equation approaches in the limit of vanishing anharmonicity. To speed up the calculation, we will use the Tersoff potential to obtain the second- and third-order force constants. We will do this with the following simple script:

```python
import numpy as np
from quippy.potential import Potential
from ase import Atoms
import ase.io
from ase.eos import calculate_eos
from ase.units import kJ
from ase.phonons import Phonons as AsePhonons
from ase.constraints import ExpCellFilter
from ase.optimize import BFGS, QuasiNewton
import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.Structure
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax


# This function will use ASE to give us a starting
# guess for dynamical matrices
def get_starting_dynamical_matrices(structure_filename,
        potential, supercell):
    atoms = ase.io.read(structure_filename)
    atoms.set_calculator(potential)
    ecf = ExpCellFilter(atoms)
    qn = QuasiNewton(ecf)
    qn.run(fmax=0.0005)

    structure = CC.Structure.Structure()
    structure.generate_from_ase_atoms(atoms, get_masses = True)
    dyn = CC.Phonons.compute_phonons_finite_displacements(structure,
    potential, supercell = supercell)
    dyn.Symmetrize()
    dyn.ForcePositiveDefinite()

    eos = calculate_eos(atoms)
    v0, e0, B = eos.fit()
    bulk = B / kJ * 1.0e24

    return dyn, bulk


# Our input variables
temperature = 100.0
nconf = 1000
max_pop = 1000

# Load in Tersoff potential
pot = Potential('IP Tersoff', param_filename=
'../06_the_SSCHA_with_machine_learning_potentials/ip.parms.Tersoff.xml')
supercell = tuple((4*np.ones(3, dtype=int)).tolist())
dyn, bulk = get_starting_dynamical_matrices(
'../06_the_SSCHA_with_machine_learning_potentials/POSCAR', pot, supercell)


# Generate the ensemble and the minimizer objects
ensemble = sscha.Ensemble.Ensemble(dyn, T0=temperature,
supercell = dyn.GetSupercell())
ensemble.generate(N = nconf)
minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
minimizer.min_step_dyn = 0.1
minimizer.kong_liu_ratio = 0.5
minimizer.meaningful_factor = 0.001
minimizer.max_ka = 1000

# Relax structure
relax = sscha.Relax.SSCHA(minimizer, ase_calculator = pot,
N_configs = nconf, max_pop = max_pop, save_ensemble = True)
relax.vc_relax(static_bulk_modulus = bulk,
ensemble_loc = "directory_of_the_ensemble")

# Generate ensemble for third-order FC with the relaxed dynamical matrices
new_ensemble = sscha.Ensemble.Ensemble(relax.minim.dyn, T0=temperature,
supercell = relax.minim.dyn.GetSupercell())
new_ensemble.generate(N = nconf*5)
new_ensemble.compute_ensemble(pot, compute_stress = True,
stress_numerical = False, cluster = None, verbose = True)
# We minimize the free energy with this new ensemble
new_minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(new_ensemble)
new_minimizer.minim_struct = False
new_minimizer.set_minimization_step(0.1)
new_minimizer.meaningful_factor = 0.001
new_minimizer.max_ka = 10000
new_minimizer.init()
new_minimizer.run()
new_minimizer.dyn.save_qe('final_dyn')
# Update weights with a new dynamical matrice
new_ensemble.update_weights(new_minimizer.dyn, temperature)

# Calculate Hessian and the third order tensor (return_d3 = True)
dyn_hessian, d3_tensor = new_ensemble.get_free_energy_hessian(include_v4 = False,
    get_full_hessian = True, return_d3 = True)
np.save("d3.npy", d3_tensor)
dyn_hessian.save_qe('hessian_dyn')
```

Here we used a 4x4x4 supercell. **You will need to converge results with respect to the size of the supercell**. A good check for convergence could be the decay of the second- and third-order force constants with distance. Now that we have second- and third-order force constants, we can calculate the lattice thermal conductivity. For this, we provide the following script:

```python
from __future__ import print_function
from __future__ import division

import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.ThermalConductivity
import time

dyn_prefix = 'final_dyn'
nqirr = 8

SSCHA_TO_MS = cellconstructor.ThermalConductivity.SSCHA_TO_MS
RY_TO_THZ = cellconstructor.ThermalConductivity.SSCHA_TO_THZ
dyn = CC.Phonons.Phonons(dyn_prefix, nqirr)

supercell = dyn.GetSupercell()

fc3 = CC.ForceTensor.Tensor3(dyn.structure,
dyn.structure.generate_supercell(supercell), supercell)

d3 = np.load("d3.npy")
fc3.SetupFromTensor(d3)
fc3 = CC.ThermalConductivity.centering_fc3(fc3)

mesh = [10,10,10]
smear = 0.03/RY_TO_THZ

tc = CC.ThermalConductivity.ThermalConductivity(dyn, fc3,
kpoint_grid = mesh, scattering_grid = mesh, smearing_scale = None,
smearing_type = 'constant', cp_mode = 'quantum', off_diag = False)

temperatures = np.linspace(200,1200,10,dtype=float)
start_time = time.time()
tc.setup_harmonic_properties(smear)
tc.write_harmonic_properties_to_file()

tc.calculate_kappa(mode = 'SRTA', temperatures = temperatures,
write_lifetimes = True, gauss_smearing = True, offdiag_mode = 'wigner',
kappa_filename = 'Thermal_conductivity_SRTA', lf_method = 'fortran-LA')

print('Calculated SSCHA kappa in: ', time.time() - start_time)

tc.calculate_kappa(mode = 'GK', write_lineshapes=False,
ne = 1000, temperatures = temperatures,
kappa_filename = 'Thermal_conductivity_GK')

print('Calculated SSCHA kappa in: ', time.time() - start_time)
# Save ThermalConductivity object for postprocessing.
tc.save_pickle()
```

Important parts of the script are:

- We define the mesh on which we calculate phonon properties to be the same as the mesh on which we calculate scattering processes (variable `mesh`). This does not have to be the case. In most cases, `scattering_grid` can be much smaller than `kpoint_grid`. **Converge your results with respect to both grids.**
- We use a smearing approach to satisfy energy conservation laws. There are two methods: constant and adaptive. In the case of `smearing_type = 'constant'`, we have to provide the smearing value in **Ry** as the argument to the `setup_harmonic_properties` function. If we choose adaptive smearing, the smearing constant will be different for different phonon modes. We can still define a global variable `smearing_scale` with which we multiply precomputed smearing constants. `smearing_scale = 1.0` works pretty well in most cases. **Converge your results with respect to smearing variables.**
- The `off_diag` variable defines whether we are doing the calculation with what was termed *coherent transport*. This will be important for highly anharmonic materials with large bunching of phonon modes.
- The function `calculate_kappa` does most of the work. Here we describe the main options:

  - `mode` defines which method to use to calculate lattice thermal conductivity. Options are **SRTA**, which is the Boltzmann transport equation solution in the single relaxation time approximation, and **GK** ([Dangic et al.](https://www.nature.com/articles/s41524-021-00523-7)), which is the Green-Kubo method that uses phonon spectral functions instead of phonon lifetimes. These two modes should give similar results in low-anharmonicity materials, but differ in strongly anharmonic ones.
  - `gauss_smearing` defines how we treat energy conservation in the calculation of self-energy. If **True**, it will use Gaussian functions; if **False**, it will use Lorentzian functions. In the case of Gaussian smearing, the real part of the self-energy is calculated using the Kramers-Kronig transformation.
  - `offdiag_mode` defines how we calculate coherent transport if `mode = 'SRTA'`. Two options: **wigner** ([Simoncelli et al.](https://www.nature.com/articles/s41567-019-0520-x)) and **gk** ([Isaeva et al.](https://www.nature.com/articles/s41467-019-11572-4)). If `mode` is **GK**, coherent transport is included naturally.
  - `lf_method` defines how lifetimes are calculated when `mode = 'SRTA'`. In short, you want to keep *fortran-*, and then add *LA* or *P*. These should give more or less the same results. An additional option is *SC*, where we solve phonon lifetimes self-consistently, meaning we account for phonon lineshifts.
  - `ne` defines the number of frequency steps if we are calculating phonon lineshapes. This is also important in the case of `lf_method = 'SC'` because we solve the self-consistent equation on a grid of frequency values, linearly interpolating the real and imaginary parts. Larger is better. **Converge your results with respect to ne.**

This calculation should take a few minutes. The results are saved in **kappa_filename**.

> **Question:**
>
> If we check the results, we see that the *SRTA* and *GK* results are different. Why? How can we improve this calculation?

## Lattice Thermal Conductivity of GeTe

As a second example, we will calculate the lattice thermal conductivity of GeTe. GeTe is a highly anharmonic material with a phase transition from rhombohedral to cubic at around 700 K. This means its lattice thermal conductivity is very low. Additionally, it should show differences between the *SRTA* and *GK* methods.

For SSCHA minimization, we can calculate atomic properties using the [Gaussian Approximation Potential](https://archive.materialscloud.org/record/2021.42) developed for this material. However, in the interest of time, we have provided the dynamical matrices calculated at 0 K and the third-order force constants in the folder `09_Thermal_conductivity_calculations_with_the_SSCHA`.

> **Exercise:**
>
> Calculate the lattice thermal conductivity of GeTe up to 1200 K (sample temperature from 300 K every 200 K). Is there a difference between the *GK* and *SRTA* methods?

> **Exercise:**
>
> Check whether coherent transport has an influence on thermal conductivity in this material system.

Finally, if we want to do some post-processing, we can load the previously saved `ThermalConductivity` object and access all previously calculated data. For example, we can calculate the phonon density of states computed with auxiliary force constants and the one computed with phonon lineshapes:

```python
from __future__ import print_function
from __future__ import division

# Import the modules to read the dynamical matrix
import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.ThermalConductivity
import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

tc = CC.ThermalConductivity.load_thermal_conductivity()

# See at which temperatures we calculated stuff
tc.what_temperatures()

key = list(tc.lineshapes.keys()) # Get Ts for lineshapes

# DOS calculated from auxiliary force constants
harm_dos = tc.get_dos()
# Temperature dependent DOS calculated from lineshapes
# first two arrays are raw data
# second two is gaussian smoothed results \
#for the distance between energy points de
anharm_dos = tc.get_dos_from_lineshapes(float(key[-1]), de = 0.1)


# Plot results
fig = plt.figure(figsize=(6.4, 4.8))
gs1 = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs1[0, 0])
ax.plot(harm_dos[0], harm_dos[1], 'k-',
zorder=0, label = 'Harmonic')
ax.plot(anharm_dos[0], anharm_dos[1], 'r-',
zorder=0, label = 'Anharmonic raw @ ' + key[-1] + ' K')
ax.plot(anharm_dos[2], anharm_dos[3], 'b-',
zorder=0, label = 'Anharmonic smooth @ ' + key[-1] + ' K')
ax.set_xlabel('Frequency (THz)')
ax.set_ylabel('Density of states')
ax.legend(loc = 'upper right')
ax.set_ylim(bottom = 0.0)
fig.savefig('test.pdf')
plt.show()
```

Additionally, if we want to check a specific phonon lineshape (for example at the $\Gamma$ point), we can do it with a bit of hacking:

```python
for iqpt in range(tc.nkpt):
    if(np.linalg.norm(tc.k_points[iqpt]) == 0.0):
        break

energies = np.arange(len(tc.lineshapes[key[-1]][0,0]),
dtype=float)*tc.delta_omega + tc.delta_omega
tc.write_lineshape('Lineshape_at_Gamma',
tc.lineshapes[key[-1]][iqpt], iqpt, energies, 'no')
```
