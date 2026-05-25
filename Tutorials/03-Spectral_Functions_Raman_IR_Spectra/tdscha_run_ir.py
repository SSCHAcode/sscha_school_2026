import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble

import tdscha, tdscha.DynamicalLanczos

import sys, os

TEMPERATURE = 450 # K
NQIRR = 4 # Irreducible q points of the dynamical matrix

# Info about the dynamical matrix and the ensemble
ORIGINAL_DYN = "sscha_auxiliary_dyn_"
FINAL_DYN = "sscha_converged_dyn_"
POP_ID = 100

def compute_ir():
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
    final_dyn.ReadInfoFromESPRESSO("eff_charges_raman.pho")

    # Update the ensemble weights on the converged dynamical matrix
    ensemble.update_weigths(FINAL_DYN, TEMPERATURE)

    # Initialize the TD-SCHA Lanczos algorithm
    lanczos = tdscha.DynamicalLanczos.Lanczos(ensemble)
    lanczos.init()

    # Let us define which level of anharmonicity we want
    lanczos.ignore_v3 = False # Add bubble contribution if false
    lanczos.ignore_v4 = False # Add RPA resummation if false

    # If both v3 and v4 are ignored (both true), we get the 'harmonic' spectrum
    # on the sscha auxiliary frequencies

    # Define the reponse function to observe
    # In this case IR with polarization along the x axis
    polarization = np.array([1.0, 0.0, 0.0])
    lanczos.prepare_ir(pol_vec = polarization)


