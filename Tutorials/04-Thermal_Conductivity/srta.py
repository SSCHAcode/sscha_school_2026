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
numconf = 2000

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
    smearing_scale=1.0,
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
    kappa_filename="Thermal_conductivity_SRTA_" + str(numconf),
    lf_method="fortran-LA",
)
