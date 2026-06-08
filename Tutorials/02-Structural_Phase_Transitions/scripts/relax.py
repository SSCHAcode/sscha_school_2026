import os
import numpy as np
import warnings
import cellconstructor as CC
import cellconstructor.Phonons
import sscha
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
TEMPERATURE = 0

# Number of irreducible q-points
NQIRR = 3

# Number of configurations used in the SSCHA ensemble
N_CONFIGS = 126

# Prefix of the starting dynamical matrices
PATH_STARTDYN="../../Materials/tutorial_02/"
START_DYN = f"{PATH_STARTDYN}ffield_dynq_"

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


# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

output_dir = f"RELAX/T_{TEMPERATURE}"
os.makedirs(output_dir, exist_ok=True)


# ============================================================
# LOAD STARTING DYNAMICAL MATRICES
# ============================================================

start_dyn = CC.Phonons.Phonons(START_DYN, NQIRR)

# If the starting dynamical matrices are not positive definite
# (i.e. imaginary phonon frequencies are present),
# enforce positive definiteness and reimpose crystal symmetries
# and the acoustic sum rule.

w, pols = start_dyn.DiagonalizeSupercell()

w_sorted = np.sort(w)
tol = -1e-4

if np.min(w_sorted[3:]) < tol:
    print("Non-acoustic imaginary phonon modes detected")
    start_dyn.ForcePositiveDefinite()
    start_dyn.Symmetrize()


# ============================================================
# CREATE SSCHA ENSEMBLE
# ============================================================

ensemble = sscha.Ensemble.Ensemble(
    start_dyn,
    TEMPERATURE
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

relax.minim.dyn.save_qe(
    os.path.join(output_dir, "sscha_dyn_")
)
