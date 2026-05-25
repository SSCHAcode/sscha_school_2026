import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import numpy as np

import tdscha, tdscha.DynamicalLanczos

import sys, os

TEMPERATURE = 300 # K
NQIRR = 3 # Irreducible q points of the dynamical matrix

# Info about the dynamical matrix and the ensemble
ORIGINAL_DYN = "sscha_auxiliary_tetra_"
FINAL_DYN = "sscha_converged_tetra_dyn_"
POP_ID = 200

def compute_raman():
    # Load the original ensemble
    dyn_original = CC.Phonons.Phonons(ORIGINAL_DYN, NQIRR)
    ensemble = sscha.Ensemble.Ensemble(dyn_original, TEMPERATURE)
    ensemble.load_bin("data", POP_ID)

    # Lets load the final converged dynamical matrix
    final_dyn = CC.Phonons.Phonons(FINAL_DYN, NQIRR)
    
    # To prepare the Raman, we need to load the Raman tensor
    # Load them from quantum espresso ph.x output
    final_dyn.ReadInfoFromESPRESSO("dielectric_calc_tetra.pho")

    # Update the ensemble weights on the converged dynamical matrix
    ensemble.update_weights(final_dyn, TEMPERATURE)

    # Initialize the TD-SCHA Lanczos algorithm
    lanczos = tdscha.DynamicalLanczos.Lanczos(ensemble, lo_to_split = None)
    lanczos.init()

    # Let us define which level of anharmonicity we want
    lanczos.ignore_v3 = False # Add bubble contribution if false
    lanczos.ignore_v4 = True # Add RPA resummation if false (a factor 2 slower in speed - no extra memory)

    # If both v3 and v4 are ignored (both true), we get the 'harmonic' spectrum
    # on the sscha auxiliary frequencies

    # Define the reponse function to observe
    # In this case IR with polarization along the x axis
    polarization_in = np.array([1.0, 0.0, 0.0])
    polarization_out = np.array([1.0, 0.0, 0.0])
    lanczos.prepare_raman(pol_vec_in = polarization_in,
                          pol_vec_out = polarization_out
                          )

    # Run the lanczos algorithm for 40 steps
    lanczos.run_FT(50)

    # Save the final result in binary
    lanczos.save_status("Raman_xx.npz")

if __name__ == "__main__":
    # Change the working directory to the one containing this script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    compute_raman()



