from __future__ import division, print_function

import time

import numpy as np

import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Phonons
import cellconstructor.ThermalConductivity

# Use the dynamical matrix and third-order force constants obtained with
# 10000 configurations.
dyn_prefix = "final_dyn_2000_"
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
d3 = np.load("d3_2000.npy")
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
ne = 6000

tc.calculate_kappa(
    mode="GK",
    temperatures=temperatures,
    write_lifetimes=False,
    gauss_smearing=False,
    ne=ne,
    offdiag_mode="perturbative",
    kappa_filename="Thermal_conductivity_GK",
    lf_method="fortran-LA",
)

# Save the ThermalConductivity object for post-processing.
tc.save_pickle()
