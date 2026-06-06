import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import tdscha, tdscha.DynamicalLanczos

import numpy as np
import sys, os

TEMPERATURE = 450
NQIRR = 4

# Info about the dynamical matrix and the ensemble
ORIGINAL_DYN = "sscha_auxiliary_dyn_"
FINAL_DYN = "sscha_converged_dyn_"
POP_ID = 100

def compute_hessian():
    # Load the original ensemble
    dyn_original = CC.Phonons.Phonons(ORIGINAL_DYN, NQIRR)
    ensemble = sscha.Ensemble.Ensemble(dyn_original, TEMPERATURE)
    ensemble.load_bin("data", POP_ID)

    # Lets load the final converged dynamical matrix
    final_dyn = CC.Phonons.Phonons(FINAL_DYN, NQIRR)
    
    # To prepare the IR or Raman, we need 
    # IR : dielectric tensor and Born effective charges
    # Raman : Raman tensor
    # Load them from quantum espresso ph.x output
    final_dyn.ReadInfoFromESPRESSO("dielectric_calc.pho")

    # Update the ensemble weights on the converged dynamical matrix
    ensemble.update_weights(final_dyn, TEMPERATURE)

    # Now we compute the free energy Hessian using the standard way
    hessian = ensemble.get_free_energy_hessian(include_v4 = False)
    hessian.save_qe("hessian")

    # Initialize the TD-SCHA Lanczos algorithm
    lanczos = tdscha.DynamicalLanczos.Lanczos(ensemble)
    lanczos.init()

    # Let us define which level of anharmonicity we want
    lanczos.ignore_v3 = False # Add bubble contribution if false
    lanczos.ignore_v4 = True # Add RPA resummation if false (a factor 2 slower in speed - no extra memory)

    # Define the perturbation of the lowest phonon mode
    lanczos.prepare_mode(0) # 0 is the lowest auxiliary phonon in the supercell

    # Run the lanczos algorithm for 40 steps
    lanczos.run_FT(40)

    # Save the final result in binary
    lanczos.save_status("hessian.npz")

    # Plot the convergence of the Hessian with
    # tdscha-convergence-analysis hessian.npz
    # Or extract the phonon mode from the Green's Function
    # w^2 = 1 / G(0)^{-1}
    green_function = lanczos.get_green_function_continued_fraction([0])[0]
    w_squared = 1.0 / np.real(green_function)
    w_hessian = np.sqrt(np.abs(w_squared)) * np.sign(w_squared)
    print("The Hessian frequency of the lowest SSCHA auxiliary mode is")
    print(" w = {:.4f}".format(w_hessian * CC.Units.RY_TO_CM))


