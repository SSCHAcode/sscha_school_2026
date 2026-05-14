# Script to perform a simple sscha relaxation
import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax

from quippy.potential import Potential

import sys, os

TEMPERATURE = 450 # K
NQIRR = 4
START_DYN = "../../Materials/CsPbI3_cubic_harmonic_"
POTENTIAL = "../../Materials/GAP_1.xml"
N_CONFIGS = 1024


def sscha_relax(temperature = TEMPERATURE):
    # Load the harmonic dynamical matrix
    dyn = CC.Phonons.Phonons(START_DYN, NQIRR)

    # Force positive phonons
    dyn.ForcePositiveDefinite()

    # Impose symmetries and ASR
    dyn.Symmetrize()

    # Load the interatomic Potential for CsPbI3
    calc = Potential("IP GAP", param_filename=POTENTIAL)

    # Run the sscha
    ensemble = sscha.Ensemble.Ensemble(dyn, temperature)
    minimizer = sscha.SchaMinimizer.SchaMinimizer(ensemble)
    relax = sscha.Relax.Relax(minimizer, calc, N_configs = N_CONFIGS)

    relax.relax()

    # Save the final ensemble
    relax.minim.ensemble.save_bin("data", 1)



# Change the working directory to the one containing this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))


