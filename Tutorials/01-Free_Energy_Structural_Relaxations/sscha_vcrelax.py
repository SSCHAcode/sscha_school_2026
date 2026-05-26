import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax

from quippy.potential import Potential

import sys, os

import ase
import numpy as np

import warnings

# NumPy moved ComplexWarning in newer versions. This block keeps the script
# compatible with both older and newer NumPy releases.
try:
    ComplexWarning = np.exceptions.ComplexWarning
except AttributeError:
    ComplexWarning = np.ComplexWarning

warnings.filterwarnings("ignore", category=ComplexWarning)

TEMPERATURE = 450 # K
NQIRR = 4
START_DYN = "../../Materials/CsPbI3_cubic_harmonic_"
POTENTIAL = "../../Materials/GAP_1.xml"
N_CONFIGS = 50 


# Load the harmonic dynamical matrix
dyn = CC.Phonons.Phonons(START_DYN, NQIRR)

# Force positive phonons
dyn.ForcePositiveDefinite()

# Impose symmetries and ASR
dyn.Symmetrize()

# Load the interatomic Potential for CsPbI3
calc = Potential("IP GAP", param_filename=POTENTIAL)

# Run the sscha
ensemble = sscha.Ensemble.Ensemble(dyn, TEMPERATURE)
minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
minimizer.meaningful_factor = 0.000001
minimizer.set_minimization_step(0.001)
minimizer.kong_liu_ratio = 0.2
relax = sscha.Relax.SSCHA(minimizer, calc, N_configs = N_CONFIGS)

# Setup to save the minimization details every good minimization step
# Allows to plot minimization results
ioinfo = sscha.Utilities.IOInfo()
ioinfo.SetupSaving("minimization_data")
relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)

# We perform a variable cell relaxation with the target pressure in GPaf 
relax.vc_relax(target_press = 0)

# Save the final ensemble
relax.minim.ensemble.save_bin("data", 1)
relax.minim.dyn.save_qe("sscha_auxiliary_dyn_")


    


