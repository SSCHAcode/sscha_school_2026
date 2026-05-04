# Hands-on Session 1 - First SSCHA Simulations: Free Energy and Structural Relaxations

In this hands-on session, we provide ready-to-use examples to set up your first SSCHA calculation.

## The Free Energy of Gold: A Simulation in the NVT Ensemble

This simple tutorial explains how to set up an SSCHA calculation starting from just the structure, in this case a CIF file downloaded from the [Materials Project](https://materialsproject.org/materials/mp-81/) database, which we can find in the `01_First_SSCHA_simulations` directory.

Starting from the Gold structure in the primitive cell, to run the SSCHA we need:

1. Compute the harmonic phonons (dynamical matrix)
2. Remove imaginary frequencies (if any)
3. Run the SSCHA

At the very beginning, we simply import the SSCHA libraries, CellConstructor, the math libraries, and the force field. This is done in Python with the `import` statements.

```python
# Import the sscha code
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

# Import the cellconstructor library to manage phonons
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.Structure, cellconstructor.calculators

# Import the force field of Gold
import ase, ase.calculators
from ase.calculators.emt import EMT

# Import numerical and general purpose libraries
import numpy as np, matplotlib.pyplot as plt
import sys, os
```

The first thing we do is initialize a CellConstructor structure from the CIF file downloaded from the material database (`Au.cif`). We initialize the EMT calculator from ASE and relax the structure:

```python
gold_structure = CC.Structure.Structure()
gold_structure.read_generic_file("Au.cif")

# Get the force field for gold
calculator = EMT()

# Relax the gold structure
relax = CC.calculators.Relax(gold_structure, calculator)
gold_structure_relaxed = relax.static_relax()
```

In the case of Gold, the relaxation is unnecessary, as it is an FCC structure with the Fm-3m symmetry group and 1 atom per primitive cell. This means the atomic positions have no degrees of freedom, so the relaxation will end before it even starts.

Next, we perform the harmonic phonon calculation using CellConstructor and a finite displacement approach:

```python
gold_harmonic_dyn = CC.Phonons.compute_phonons_finite_displacements(gold_structure_relaxed, calculator, supercell = (4,4,4))

# Impose the symmetries and
# save the dynamical matrix in the quantum espresso format
gold_harmonic_dyn.Symmetrize()
gold_harmonic_dyn.save_qe("gold_harmonic_dyn")
```

The method `compute_phonons_finite_displacements` is documented in the CellConstructor guide. It requires the structure (in this case `gold_structure_relaxed`), the force field (`calculator`), and the supercell for the calculation. Here we use a 4x4x4 supercell (equivalent to 64 atoms). This may not be sufficient to converge all properties, especially at very high temperature, but it is a good starting point.

Note that `compute_phonons_finite_displacements` works in parallel with MPI. Therefore, if the script is executed with `mpirun -np 16 python myscript.py`, it will split the calculations of the finite displacements across 16 processors. You need to have mpi4py installed. However, in this case, due to the high symmetry, only one calculation is required to obtain the harmonic dynamical matrix, so parallelization is unnecessary.

After computing the harmonic phonons in `gold_harmonic_dyn`, we impose the correct symmetrization and the acoustic sum rule with the `Symmetrize` method, and save the result in the Quantum ESPRESSO format with `save_qe`.

We are ready to submit the SSCHA calculation in the NVT ensemble.

The important parameters are:

- The temperature
- The number of random configurations in the ensemble
- The maximum number of iterations

These parameters are almost self-explanatory. In contrast with Molecular Dynamics (MD) or Metropolis-Monte Carlo (MC) calculations, where the equilibrium probability distribution is sampled through a dynamical evolution of the structure, the SSCHA encodes the whole probability distribution as an analytical function. Therefore, to compute properties, we can generate on the fly the configurations that sample the equilibrium distribution.

![Workflow of the SSCHA objects for a free energy minimization.](figures_01/Diagram1.png)

The code that sets up and performs the SSCHA is the following:

```python
TEMPERATURE = 300
N_CONFIGS = 50
MAX_ITERATIONS = 20

# Initialize the random ionic ensemble
ensemble = sscha.Ensemble.Ensemble(gold_harmonic_dyn, TEMPERATURE)

# Initialize the free energy minimizer
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
minim.set_minimization_step(0.01)

# Initialize the NVT simulation
relax = sscha.Relax.SSCHA(minim, calculator, N_configs = N_CONFIGS,
max_pop = MAX_ITERATIONS)

# Define the I/O operations
# To save info about the free energy minimization after each step
ioinfo = sscha.Utilities.IOInfo()
ioinfo.SetupSaving("minim_info")
relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)

# Run the NVT simulation
relax.relax(get_stress = True)

# Save the final dynamical matrix
relax.minim.dyn.save_qe("sscha_T300_dyn")
```

In the previous code we defined the main objects needed to run the simulation:

- `ensemble` (`sscha.Ensemble.Ensemble`) represents the ensemble of ionic configurations. We initialize it with the dynamical matrix (which determines how much atoms fluctuate around the centroids) and the temperature.
- `minim` (`sscha.SchaMinimizer.SSCHA_Minimizer`) performs the free energy minimization. It contains all the information regarding the minimization algorithm, such as the initial timestep (set to 0.01 here). You can avoid setting the timestep, as the code will automatically guess the best value.
- `relax` (`sscha.Relax.SSCHA`) automates the generation of ensembles, calculation of energies and forces, and the free energy minimization to perform an NVT or NPT calculation. To initialize it, we pass the `minim` (which contains the ensemble with the temperature), the force field (`calculator`), the number of configurations `N_configs`, and the maximum number of iterations.

In this example, most of the time is spent in the minimization. However, if we replace the force field with ab initio DFT, the time to run the minimization is negligible compared to the time required to compute energies and forces on the ensemble configurations. The total (maximum) number of energy/forces calculations is equal to the number of configurations times the number of iterations (passed through the `max_pop` argument).

The calculation is submitted with `relax.relax()`. However, before running the calculation we introduce another object, the `IOInfo`. This tells `relax` to save information about the free energy, its gradient, and the anharmonic phonon frequencies during the minimization in the files `minim_info.dat` and `minim_info.freqs`. It is not mandatory to use them, but it is very useful as it allows visualization of the minimization while it is running.

The full input file is:

```python
# Import the sscha code
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax, sscha.Utilities

# Import the cellconstructor library to manage phonons
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.Structure, cellconstructor.calculators

# Import the force field of Gold
import ase, ase.calculators
from ase.calculators.emt import EMT

# Import numerical and general purpose libraries
import numpy as np, matplotlib.pyplot as plt
import sys, os


"""
Here we load the primitive cell of Gold from a cif file.
And we use CellConstructor to compute phonons from finite differences.
The phonons are computed on a q-mesh 4x4x4
"""

gold_structure = CC.Structure.Structure()
gold_structure.read_generic_file("Au.cif")

# Get the force field for gold
calculator = EMT()

# Relax the gold structure (useless since for symmetries it is already relaxed)
relax = CC.calculators.Relax(gold_structure, calculator)
gold_structure_relaxed = relax.static_relax()

# Compute the harmonic phonons
# NOTE: if the code is run with mpirun, the calculation goes in parallel
gold_harmonic_dyn = CC.Phonons.compute_phonons_finite_displacements(gold_structure_relaxed, calculator, supercell = (4,4,4))

# Impose the symmetries and
# save the dynamical matrix in the quantum espresso format
gold_harmonic_dyn.Symmetrize()
gold_harmonic_dyn.save_qe("harmonic_dyn")


# If the dynamical matrix has imaginary frequencies, remove them
gold_harmonic_dyn.ForcePositiveDefinite()

"""
gold_harmonic_dyn is ready to start the SSCHA calculation.

Now let us initialize the ensemble, and the calculation at 300 K.
We will run a NVT calculation, using 100 configurations at each step
"""

TEMPERATURE = 300
N_CONFIGS = 50
MAX_ITERATIONS = 20

# Initialize the random ionic ensemble
ensemble = sscha.Ensemble.Ensemble(gold_harmonic_dyn, TEMPERATURE)

# Initialize the free energy minimizer
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
minim.set_minimization_step(0.01)

# Initialize the NVT simulation
relax = sscha.Relax.SSCHA(minim, calculator, N_configs = N_CONFIGS,
max_pop = MAX_ITERATIONS)

# Define the I/O operations
# To save info about the free energy minimization after each step
ioinfo = sscha.Utilities.IOInfo()
ioinfo.SetupSaving("minim_info")
relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)


# Run the NVT simulation (save the stress to compute the pressure)
relax.relax(get_stress = True)

# Save the final dynamical matrix
# And print in stdout the info about the minimization
relax.minim.finalize()
relax.minim.dyn.save_qe("sscha_T{}_dyn".format(TEMPERATURE))
```

Now save the file as `sscha_gold.py` and execute it with:

```bash
$ python sscha_gold.py > output.log
```

And that is it. The code will probably take a few minutes on a standard laptop computer.

**Congratulations!** You ran your first SSCHA simulation!

If you open a new terminal in the same directory as the SSCHA submission, you can plot the information during the minimization. Starting from version 1.2, we provide visualization utilities installed together with the SSCHA. Simply type:

```bash
$ sscha-plot-data.py minim_info
```

You will see two windows.

![Minimization data of Gold.](figures_01/gold_minim.png)

The figure above shows all the minimization data. In the top-left panel, we see the free energy. As expected, it decreases (since the SSCHA is minimizing it). You can see that at certain values of the steps there are discontinuities. These occur when the code realizes that the ensemble on which it is computing is no longer good and a new one is generated. The quality of an ensemble is determined by the Kong-Liu effective sample size (bottom-left). When it reaches 0.5 of its initial value (equal to the number of configurations), the ensemble is extracted again and a new iteration starts. You can see that in the last iteration, the code stops before getting to 25 ($0.5 \cdot 50$). This means that the code converged properly: the gradient reached zero while the ensemble was still good.

On the right side you see the free energy gradients, which must go to zero to converge. The top-right is the gradient of the SSCHA dynamical matrix, while on bottom-right there is the gradient of the average atomic positions. Indeed, since the gold atomic positions are all fixed by symmetries, it is always zero (but it will be different from zero in more complex systems).

![All the SSCHA phonon frequencies as a function of the step in the NVT simulation.](figures_01/frequency_minim.png)

The figure above represents the evolution of the SSCHA phonon frequencies. Here, all the frequencies in the supercell (at each q-point commensurate with the calculation) are shown.

> **Note:** The SSCHA auxiliary frequencies in the figure above are not the real frequencies observed in experiments, but are rather linked to the average displacements of atoms along that mode.

By looking at how they change, you can get an idea of which phonon modes are more affected by anharmonicity. In this case, it is evident that Gold is strongly anharmonic and that the temperature makes all the phonon frequencies harder.

At the end of the simulation, the code writes the final dynamical matrix in the Quantum ESPRESSO file format: `sscha_T300_dynX` where X goes over the number of irreducible q-points.

In the next section, we analyze in detail each section of the script to provide a bit more insight into the simulation, and a guide to modify it to fit your needs and submit your own system.

## Plot the Phonon Dispersion

Now that the SSCHA minimization has ended, we can compare the harmonic and anharmonic phonon dispersion of Gold.

For this purpose, we can simply run a script like the following. You can find a copy of this script already in `Examples/ThermodynamicsOfGold/plot_dispersion.py`. You can use it even in your own simulation; simply edit the value of the uppercase keyword at the beginning of the script to match your needs.

```python
# Import the CellConstructor library to plot the dispersion
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.ForceTensor

# Import the numerical libraries and those for plotting
import numpy as np
import matplotlib.pyplot as plt

import sys, os

# Let us define the PATH in the brilluin zone and the total number of points
PATH = "GXWXKGL"
N_POINTS = 1000

# Here we define the position of the special points
SPECIAL_POINTS = {"G": [0,0,0],
    "X": [0, .5, .5],
    "L": [.5, .5, .5],
    "W": [.25, .75, .5],
    "K": [3/8., 3/4., 3/8.]}

# The two dynamical matrix to be compared
HARM_DYN = 'harmonic_dyn'
SSCHA_DYN = 'sscha_T300_dyn'

# The number of irreducible q points
# i.e., the number of files in which the phonons are stored
NQIRR = 13


# --------------------- THE SCRIPT FOLLOWS ---------------------

# Load the harmonic and sscha phonons
harmonic_dyn = CC.Phonons.Phonons(HARM_DYN, NQIRR)
sscha_dyn = CC.Phonons.Phonons(SSCHA_DYN, NQIRR)

# Get the band path
qpath, data = CC.Methods.get_bandpath(harmonic_dyn.structure.unit_cell,
    PATH,
    SPECIAL_POINTS,
    N_POINTS)
xaxis, xticks, xlabels = data # Info to plot correclty the x axis

# Get the phonon dispersion along the path
harmonic_dispersion = CC.ForceTensor.get_phonons_in_qpath(harmonic_dyn, qpath)
sscha_dispersion = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn, qpath)

nmodes = harmonic_dyn.structure.N_atoms * 3

# Plot the two dispersions
plt.figure(dpi = 150)
ax = plt.gca()

for i in range(nmodes):
    lbl=None
    lblsscha = None
    if i == 0:
        lbl = 'Harmonic'
        lblsscha = 'SSCHA'

    ax.plot(xaxis, harmonic_dispersion[:,i], color='k',
            ls='dashed', label=lbl)
    ax.plot(xaxis, sscha_dispersion[:,i],
            color='r', label=lblsscha)

# Plot vertical lines for each high symmetry points
for x in xticks:
   ax.axvline(x, 0, 1, color = "k", lw = 0.4)
ax.axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)

# Set the x labels to the high symmetry points
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)

ax.set_xlabel("Q path")
ax.set_ylabel("Phonons [cm-1]")

ax.legend()

plt.tight_layout()
plt.savefig("dispersion.png")
plt.show()
```

If we save the script as `plot_dispersion.py` in the same directory as the calculation, we can run it with:

```bash
$ python plot_dispersion.py
```

![Comparison between the SSCHA and the harmonic phonon dispersion of Gold.](figures_01/gold_dispersion.png)

The script will plot the figure of the phonon dispersion shown above. It is quite different from the experimental one because of the poor accuracy of the force field; however, the SSCHA result is much closer to the experimental value.

> **Exercise:**
>
> Try to perform the simulation of Gold at a different temperature, then plot the SSCHA phonon dispersion as a function of temperature.
>
> How do the phonon bands behave if the temperature is increased? Do they become more rigid (energy increases) or softer?

## Running in the NPT Ensemble: Simulating Thermal Expansion

Now that you have some experience with the NVT simulation, we are ready for the next step: NPT, or relaxing the lattice.

With python-sscha it is very easy to run NPT simulation. You simply have to replace the line of the NVT script with the target pressure for the simulation:

```python
# Replace the line
# relax.relax(get_stress = True)
# with
relax.vc_relax(target_press = 0)
```

And that is all! The target pressure is expressed in GPa; in this case, 0 is ambient conditions (1 atm = 0.0001 GPa).

You can also perform NVT simulation with variable lattice parameters. In this case the system will constrain the total volume to remain constant, but the lattice parameters will be optimized (if the system is not cubic and has some degrees of freedom, which is not the case for Gold).

The NVT ensemble with variable lattice parameters (cell shape) is:

```python
# Replace the line
#    relax.vc_relax(target_press = 0)
# with
relax.vc_relax(fix_volume = True)
```

Indeed, this is an NVT simulation, so there is no need to specify the target pressure.

In the following script, we run the NPT ensemble at various temperatures, each time starting from the previous ensemble, to follow the volume thermal expansion of gold. This script assumes you have already performed the NVT calculation, so we can start from that result and avoid the harmonic calculation (it is always a good practice to start with an NVT simulation and then run NPT from the final result).

```python
# Import the sscha code
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sscha.Utilities

# Import the cellconstructor library to manage phonons
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.Structure, cellconstructor.calculators

# Import the force field of Gold
import ase, ase.calculators
from ase.calculators.emt import EMT

# Import numerical and general purpose libraries
import numpy as np, matplotlib.pyplot as plt
import sys, os


# Define the temperature range (in K)
T_START = 300
T_END = 1000
DT = 50

N_CONFIGS = 50
MAX_ITERATIONS = 10

# Import the gold force field
calculator = EMT()

# Import the starting dynamical matrix (final result of get_gold_free_energy.py)
dyn = CC.Phonons.Phonons("sscha_T300_dyn", nqirr = 13)

# Create the directory on which to store the output
DIRECTORY = "thermal_expansion"
if not os.path.exists(DIRECTORY):
    os.makedirs("thermal_expansion")

# We cycle over several temperatures
t = T_START


volumes = []
temperatures = []
while t <= T_END:
    # Change the temperature
    ensemble = sscha.Ensemble.Ensemble(dyn, t)
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    minim.set_minimization_step(0.1)

    relax = sscha.Relax.SSCHA(minim, calculator, N_configs = N_CONFIGS,
                              max_pop = MAX_ITERATIONS)

    # Setup the I/O
    ioinfo = sscha.Utilities.IOInfo()
    ioinfo.SetupSaving( os.path.join(DIRECTORY, "minim_t{}".format(t)))
    relax.setup_custom_functions( custom_function_post = ioinfo.CFP_SaveAll)


    # Run the NPT simulation
    relax.vc_relax(target_press = 0)

    # Save the volume and temperature
    volumes.append(relax.minim.dyn.structure.get_volume())
    temperatures.append(t)

    # Start the next simulation from the converged value at this temperature
    relax.minim.dyn.save_qe(os.path.join(DIRECTORY,
                            "sscha_T{}_dyn".format(t)))
    dyn = relax.minim.dyn

    # Print in standard output
    relax.minim.finalize()

    # Update the temperature
    t += DT

    # Save the thermal expansion
    np.savetxt(os.path.join(DIRECTORY, "thermal_expansion.dat"),
               np.transpose([temperatures, volumes]),
           header = "Temperature [K]; Volume [A^3]")
```

You can run the script as always with:

```bash
$ python thermal_expansion.py
```

And... done!

This calculation is going to require a bit more time, as we run multiple SSCHA calculations at several temperatures. After it finishes, you can plot the results written in the file `thermal_expansion/thermal_expansion.dat`.

A simple script to plot the thermal expansion (and fit the volumetric thermal expansion coefficient) is the following:

```python
import numpy as np
import matplotlib.pyplot as plt

import scipy, scipy.optimize


# Load all the dynamical matrices and compute volume
DIRECTORY = "thermal_expansion"
FILE = os.path.join(DIRECTORY, "thermal_expansion.dat")

# Load the data from the final data file
temperatures, volumes = np.loadtxt(FILE, unpack = True)


# Prepare the figure and plot the V(T) from the sscha data
plt.figure(dpi = 150)
plt.scatter(temperatures, volumes, label = "SSCHA data")

# Fit the data to estimate the volumetric thermal expansion coefficient
def parabola(x, a, b, c):
    return a + b*x + c*x**2
def diff_parab(x, a, b, c):
    return b + 2*c*x

popt, pcov = scipy.optimize.curve_fit(parabola, temperatures, volumes,
                     p0 = [0,0,0])

# Evaluate the volume thermal expansion
vol_thermal_expansion = diff_parab(300, *popt) / parabola(300, *popt)
plt.text(0.6, 0.2, r"$\alpha_v = "+"{:.1f}".format(vol_thermal_expansion*1e6)+r"\times 10^6 $ K$^{-1}$",
    transform = plt.gca().transAxes)


# Plot the fit
_t_ = np.linspace(np.min(temperatures), np.max(temperatures), 1000)
plt.plot(_t_, parabola(_t_, *popt), label = "Fit")

# Adjust the plot adding labels, legend, and saving in eps
plt.xlabel("Temperature [K]")
plt.ylabel(r"Volume [$\AA^3$]")
plt.legend()
plt.tight_layout()
plt.savefig("thermal_expansion.png")
plt.show()
```

![Thermal expansion of Gold.](figures_01/thermal_expansion.png)

We report the final thermal expansion in the figure above. The volumetric expansion coefficient $\alpha_v$ is obtained from the fit thanks to the thermodynamic relation:

$$
\alpha_v = \frac{1}{V} \left(\frac{dV}{dT}\right)_P
$$

Also in this case, the result is quite off from experiments, due to the not completely realistic force field employed. To get a more realistic result, you should use ab initio calculations or a more refined force field.

## Ab Initio Calculation with the SSCHA Code

The SSCHA code is compatible with the Atomic Simulation Environment (ASE), which we used in the previous tutorial to get a fast force field for Gold.

However, ASE already provides an interface with most codes to run ab initio simulations. The simplest way of interfacing the SSCHA to another ab initio code is to directly use ASE.

The only difference is in the definition of the calculator. In the first example of this chapter, the Gold force field was defined as:

```python
import ase
from ase.calculators.emt import EMT
calculator = EMT()
```

We simply need to replace these lines with our favorite DFT code. In this example we are going to use Quantum ESPRESSO, but the procedure for VASP, CASTEP, CRYSTAL, ABINIT, SIESTA, or your favorite one is exactly the same (refer to the official documentation of ASE for instructions on how to initialize these calculators).

In the case of DFT, unfortunately, we cannot simply create the calculator in one line, like we did for the EMT force field, as we also need to provide many parameters, such as pseudopotentials, the choice of exchange-correlation functional, the cutoff of the basis set, and the k-mesh grid for Brillouin zone sampling.

In the following example, we initialize the Quantum ESPRESSO calculator for Gold.

```python
import cellconstructor.calculators
from ase.calculators.espresso import Espresso

# Initialize the DFT (Quantum Espresso) calculator for gold
# The input data is a dictionary that encodes the pw.x input file namelist
input_data = {
    'control' : {
       # Avoid writing wavefunctions on the disk
       'disk_io' : 'None',
       # Where to find the pseudopotential
       'pseudo_dir' : '.',
       'tprnfor' : True,   # Print the forces in output
       'tstress' : True    # Print the stress tensor
    },
    'system' : {
       # Specify the basis set cutoffs
       'ecutwfc' : 45,   # Cutoff for wavefunction
       'ecutrho' : 45*4, # Cutoff for the density
       # Information about smearing (it is a metal)
       'occupations' : 'smearing',
       'smearing' : 'mv',
       'degauss' : 0.03
    },
    'electrons' : {
        'conv_thr' : 1e-8
    }
}

# the pseudopotential for each chemical element
# In this case just Gold
pseudopotentials = {'Au' : 'Au_ONCV_PBE-1.0.oncvpsp.upf'}

# the kpoints mesh and the offset
kpts = (1,1,1)
koffset = (1,1,1)


# Prepare the quantum espresso calculator
calculator = Espresso(input_data = input_data, pseudopotentials = pseudopotentials,
                         kpts = kpts, koffset = koffset)
```

If you are familiar with the Quantum ESPRESSO input files, you should recognize all the options inside the `input_data` dictionary. For more options and more information, refer to the [Quantum ESPRESSO pw.x input guide](https://www.quantum-espresso.org/Doc/INPUT_PW.html).

Remember, the parameters set here are just for demonstration; remember to run appropriate convergence checks for the k-mesh, smearing, and basis set cutoffs before running the SSCHA code. Keep also in mind that this input file refers to the supercell, and the `kpts` variable should be properly rescaled if the supercell size is increased.

All the rest of the code remains the same (but here we do not compute harmonic phonons, which can be done more efficiently within Quantum ESPRESSO). Instead, we take the result obtained with EMT in the previous sections and try to relax the free energy with a fully ab initio approach.

The complete code is inside `Examples/sscha_and_dft/nvt_local.py`.

```python
# Import the sscha code
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax, sscha.Utilities

# Import the cellconstructor library to manage phonons
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.Structure, cellconstructor.calculators

# Import the DFT calculator
import cellconstructor.calculators
from ase.calculators.espresso import Espresso

# Import numerical and general purpose libraries
import numpy as np, matplotlib.pyplot as plt
import sys, os


# Initialize the DFT (Quantum Espresso) calculator for gold
# The input data is a dictionary that encodes the pw.x input file namelist
input_data = {
    'control' : {
       # Avoid writing wavefunctions on the disk
       'disk_io' : 'None',
       # Where to find the pseudopotential
       'pseudo_dir' : '.'
    },
    'system' : {
       # Specify the basis set cutoffs
       'ecutwfc' : 45,   # Cutoff for wavefunction
       'ecutrho' : 45*4, # Cutoff for the density
       # Information about smearing (it is a metal)
       'occupations' : 'smearing',
       'smearing' : 'mv',
       'degauss' : 0.03
    },
    'electrons' : {
       'conv_thr' : 1e-8
    }
}

# the pseudopotential for each chemical element
# In this case just Gold
pseudopotentials = {'Au' : 'Au_ONCV_PBE-1.0.oncvpsp.upf'}

# the kpoints mesh and the offset
kpts = (1,1,1)
koffset = (1,1,1)

# Specify the command to call quantum espresso
command = 'pw.x -i PREFIX.pwi > PREFIX.pwo'


# Prepare the quantum espresso calculator
 calculator = Espresso(input_data =  input_data,
                    pseudopotentials = pseudopotentials,
                    command = command,
                    kpts = kpts,
                    koffset = koffset)


TEMPERATURE = 300
N_CONFIGS = 50
MAX_ITERATIONS = 20
START_DYN = 'harmonic_dyn'
NQIRR = 13

# Let us load the starting dynamical matrix
gold_dyn = CC.Phonons.Phonons(START_DYN, NQIRR)

# Initialize the random ionic ensemble
ensemble = sscha.Ensemble.Ensemble(gold_dyn, TEMPERATURE)

# Initialize the free energy minimizer
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
minim.set_minimization_step(0.01)

# Initialize the NVT simulation
relax = sscha.Relax.SSCHA(minim, calculator, N_configs = N_CONFIGS,
             max_pop = MAX_ITERATIONS)

# Define the I/O operations
# To save info about the free energy minimization after each step
ioinfo = sscha.Utilities.IOInfo()
ioinfo.SetupSaving("minim_info")
relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)


# Run the NVT simulation (save the stress to compute the pressure)
relax.relax(get_stress = True)

# If instead you want to run a NPT simulation, use
# The target pressure is given in GPa.
#relax.vc_relax(target_press = 0)

# You can also run a mixed simulation (NVT) but with variable lattice parameters
#relax.vc_relax(fix_volume = True)

# Now we can save the final dynamical matrix
# And print in stdout the info about the minimization
relax.minim.finalize()
relax.minim.dyn.save_qe("sscha_T{}_dyn".format(TEMPERATURE))
```

Now you can run the SSCHA with an ab initio code! However, your calculation will probably take forever. Let us discuss parallelization and how to exploit modern HPC infrastructures to speed things up.

## Parallelization

If you actually tried to run the code from the previous section on a laptop, it would take forever. The reason is that DFT calculations are much more expensive than the SSCHA minimization. While the SSCHA minimizes the number of ab initio calculations (especially when compared with MD or PIMD), they remain the bottleneck of the computational time.

For this reason, we need an appropriate parallelization strategy to reduce the total time to run an SSCHA.

The simplest way is to call the previous Python script with MPI:

```bash
$ mpirun -np 50 python nvt_local.py > output.log
```

The code will split the configurations in each ensemble across different MPI processes. In this case, we have 50 configurations per ensemble; by splitting them into 50 processors, we run the full ensemble in parallel.

However, even a single DFT calculation on 1 processor is going to take hours, and in some cases it may even take days. Luckily, Quantum ESPRESSO (and many other software) also has internal parallelization to work with. For example, we can tell Quantum ESPRESSO to run itself in parallel on 8 processors. For this purpose, we simply need to modify the command used to run Quantum ESPRESSO in the previous script.

```python
# Lets replace
# command = 'pw.x -i PREFIX.pwi > PREFIX.pwo'
# with
command = 'mpirun -np 8 pw.x -npool 1 -i PREFIX.pwi > PREFIX.pwo'

# The command string is passed to the espresso calculator
calculator = CC.calculators.Espresso(input_data,
                    pseudopotentials,
                    command = command,
                    kpts = kpts,
                    koffset = koffset)
```

In this way, our calculations will run on 400 processors (50 processors split the ensemble times 8 processors per calculation). This is achieved by nesting MPI calls. However, only the CellConstructor calculators can nest MPI calls without raising errors. This is the reason why we imported the Espresso class from CellConstructor and not from ASE. If you want to use ASE for your calculator, you can only use the internal parallelization of the calculator by modifying the command, as ASE itself implements an MPI parallelization on I/O operations that conflicts with the python-sscha parallelization. This limitation only applies to `FileIOCalculators` from ASE (so the EMT force field is not affected and can be safely used with python-sscha parallelization).

With this setup, the full code is parallelized over 400 processors. However, the SSCHA minimization algorithm is serial, and all the time spent in the actual SSCHA minimization is wasting the large number of allocated resources. Moreover, the SSCHA code needs to be configured and correctly installed on the cluster, which may be a difficult operation due to the hybrid Fortran/Python structure.
