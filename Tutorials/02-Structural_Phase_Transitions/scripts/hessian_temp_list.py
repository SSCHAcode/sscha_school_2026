import os
import numpy as np
import warnings

import cellconstructor as CC
import cellconstructor.Phonons

import sscha.Ensemble
import sscha.SchaMinimizer
import sscha.Relax

import fforces as ff
import fforces.Calculator


# ============================================================
# Ignore NumPy ComplexWarning generated during SSCHA updates
# ============================================================

try:
    ComplexWarning = np.exceptions.ComplexWarning
except AttributeError:
    ComplexWarning = np.ComplexWarning

warnings.filterwarnings("ignore", category=ComplexWarning)


# ============================================================
# INPUT PARAMETERS
# ============================================================

# Temperature in Kelvin
TEMPERATURE = [0, 50, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300 ]

NQIRR = 3

N_CONFIGS = 10000

# Dynamical matrices used to define the harmonic reference of the toy model
PATH_MODELDYN="../../Materials/tutorial_02/"
MODEL_DYN =f"{PATH_MODELDYN}ffield_dynq_"


# ============================================================
# DEFINE THE TOY FORCE FIELD
# ============================================================

# Load the harmonic dynamical matrices
ff_dyn = CC.Phonons.Phonons(MODEL_DYN, NQIRR)

# Create the toy-model calculator
ff_calculator = ff.Calculator.ToyModelCalculator(ff_dyn)

# Select the type of anharmonic toy potential
ff_calculator.type_cal = "pbtex"

# Anharmonic parameters
ff_calculator.p3 = 0.036475
ff_calculator.p4 = -0.022
ff_calculator.p4x = -0.014


for T in TEMPERATURE:
    
    print(f"\n===== T = {T} K =====")

    # Output directory for the free - energy Hessian dynamical matrices
    
    OUTPUT_DIR = f"HESSIAN/T_{T}"
    os. makedirs (OUTPUT_DIR , exist_ok =True)

    # ============================================================
    # LOAD RELAXED SSCHA DYNAMICAL MATRICES
    # ============================================================

    # Already relaxed SSCHA dynamical matrices
    SSCHA_DYN = f"./RELAX/T_{T}/sscha_dyn_"
    sscha_dyn = CC.Phonons.Phonons(SSCHA_DYN, NQIRR)
    sscha_dyn.Symmetrize()

    # ============================================================
    # CREATE SSCHA ENSEMBLE
    # ============================================================

    ensemble = sscha.Ensemble.Ensemble(
        sscha_dyn,
        T
    )

    # ============================================================
    # SETUP THE SSCHA MINIMIZER
    # ============================================================

    minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

    # Control minimization stability and convergence
    minimizer.meaningful_factor = 0.01
    minimizer.set_minimization_step(0.001)


    # ============================================================
    # SETUP SSCHA RELAXATION
    # ============================================================

    relax = sscha.Relax.SSCHA(
        minimizer,
        ase_calculator=ff_calculator,
        N_configs=N_CONFIGS
    )


    # ============================================================
    # RUN THE SSCHA MINIMIZATION
    # ============================================================

    relax.relax()

    # ============================================================
    # SAVE FINAL SSCHA DYNAMICAL MATRICES
    # ============================================================

    relax .minim.dyn. save_qe (
        os.path.join( OUTPUT_DIR , "refined_sscha_dyn_")
        )  


    # Reweight the ensemble using the refined dynamical matrix
    ensemble.update_weights(minimizer.dyn, T)

    # Compute the free-energy Hessian
    dyn_hessian = ensemble.get_free_energy_hessian(
        include_v4=False,
        get_full_hessian=True,
        return_d3=False,
    )

    # Save the Hessian dynamical matrix
    dyn_hessian.save_qe(
        os.path.join( OUTPUT_DIR , "hessian_dyn_")
        )


    print("\nDone.")
