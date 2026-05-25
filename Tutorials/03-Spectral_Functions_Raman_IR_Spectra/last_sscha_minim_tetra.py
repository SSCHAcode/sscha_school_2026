import sys, os
import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer

from quippy.potential import Potential


TEMPERATURE = 300 # K
NQIRR = 3
START_DYN = "sscha_auxiliary_tetra_"
POTENTIAL = "../../Materials/GAP_1.xml"
N_CONFIGS = 1024
POP_ID = 200


def last_sscha_relax(temperature = TEMPERATURE):
    # Load the sscha dynamical matrix
    dyn = CC.Phonons.Phonons(START_DYN, NQIRR)

    # Load the interatomic Potential for CsPbI3
    calc = Potential("IP GAP", param_filename=POTENTIAL)

    # Generate the last esemble
    ensemble = sscha.Ensemble.Ensemble(dyn, temperature)
    ensemble.generate(N_CONFIGS)

    # Compute the energies and forces of the ensemble with the GAP potential
    ensemble.compute_ensemble(calc)

    # Compute a full sscha minimization on the new bigger ensemble
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    minim.set_minimization_step(0.02)
    minim.meaningful_factor = 0.01
    minim.run()

    # Save the final dynamical matrix and ensemble for further calculations
    minim.ensemble.save_bin("data", POP_ID)
    minim.dyn.save_qe("sscha_converged_tetra_dyn_")


if __name__ == "__main__":
    # Change the working directory to the one containing this script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Run the final sscha relaxation
    last_sscha_relax()




