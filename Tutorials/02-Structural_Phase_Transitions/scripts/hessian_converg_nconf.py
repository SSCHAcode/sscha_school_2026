import os
import numpy as np
import warnings

import cellconstructor as CC
import cellconstructor.Phonons

import sscha.Ensemble
import sscha.SchaMinimizer

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

TEMPERATURE = 0
NQIRR = 3

N_CONFIGS_LIST = [126, 250, 500, 1000, 2000, 4000, 8000]

# Already relaxed SSCHA dynamical matrices
SSCHA_DYN = "./RELAX/T_0/sscha_dyn_"

# Dynamical matrices used to define the harmonic reference of the toy model
PATH_MODELDYN="../../Materials/tutorial_02/"
MODEL_DYN =f"{PATH_MODELDYN}ffield_dynq_"

OUTPUT_DIR = f"HESSIAN_CONVERGENCE_T_{TEMPERATURE}"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# DEFINE THE TOY FORCE FIELD
# ============================================================

ff_dyn = CC.Phonons.Phonons(MODEL_DYN, NQIRR)

calculator = ff.Calculator.ToyModelCalculator(ff_dyn)
calculator.type_cal = "pbtex"
calculator.p3 = 0.036475
calculator.p4 = -0.022
calculator.p4x = -0.014


# ============================================================
# LOAD RELAXED SSCHA DYNAMICAL MATRICES
# ============================================================

sscha_dyn = CC.Phonons.Phonons(SSCHA_DYN, NQIRR)
sscha_dyn.Symmetrize()


# ============================================================
# LOOP OVER ENSEMBLE SIZE
# ============================================================

results = []

for nconf in N_CONFIGS_LIST:

    print(f"\nRunning calculation with {nconf} configurations")

    # Generate a new ensemble using the relaxed SSCHA dynamical matrix
    ensemble = sscha.Ensemble.Ensemble(
        sscha_dyn,
        T0=TEMPERATURE,
        supercell=sscha_dyn.GetSupercell(),
    )

    ensemble.generate(N=nconf)

    # Compute energies and forces for all configurations
    ensemble.compute_ensemble(
        calculator,
        compute_stress=False,
        cluster=None,
        verbose=True,
    )

    # Minimize the auxiliary dynamical matrix at fixed structure/cell
    minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

    minimizer.minim_struct = False
    minimizer.set_minimization_step(0.001)
    minimizer.meaningful_factor = 0.001
    minimizer.max_ka = 10000

    minimizer.init()
    minimizer.run()

    # Save the refined auxiliary dynamical matrix
    refined_prefix = os.path.join(OUTPUT_DIR, f"refined_dyn_{nconf}_")
    minimizer.dyn.save_qe(refined_prefix)

    # Reweight the ensemble using the refined dynamical matrix
    ensemble.update_weights(minimizer.dyn, TEMPERATURE)

    # Compute the free-energy Hessian
    dyn_hessian = ensemble.get_free_energy_hessian(
        include_v4=False,
        get_full_hessian=True,
        return_d3=False,
    )

    # Save the Hessian dynamical matrix
    hessian_prefix = os.path.join(OUTPUT_DIR, f"hessian_dyn_{nconf}_")
    dyn_hessian.save_qe(hessian_prefix)

    # Extract Gamma frequencies from the free-energy Hessian
    w_gamma, pols = dyn_hessian.DiagonalizeSupercell()
    
    # Convert frequencies to cm^-1
    w_gamma = np.real(w_gamma) * CC.Units.RY_TO_CM
    
    w_sorted = np.sort(w_gamma)
    
    print("Gamma frequencies [cm^-1]:")
    print(w_sorted)
    
    # Remove the three acoustic modes closest to zero
    idx_sorted_by_abs = np.argsort(np.abs(w_sorted))
    acoustic_indices = idx_sorted_by_abs[:3]
    
    optical_modes = np.delete(w_sorted, acoustic_indices)
    optical_modes = np.sort(optical_modes)

    lowest_optical_gamma = optical_modes[0]

    print(f"Lowest optical Gamma frequency = {lowest_optical_gamma:.6f} cm^-1")

    results.append([nconf, lowest_optical_gamma])

# ============================================================
# SAVE CONVERGENCE DATA
# ============================================================

results = np.array(results)

# Create output table:
# column 1 -> N_configs
# column 2 -> 1/sqrt(N_configs)
# column 3 -> lowest optical Gamma frequency

output_data = np.column_stack([
    results[:, 0],
    1 / np.sqrt(results[:, 0]),
    results[:, 1],
])

np.savetxt(
    os.path.join(OUTPUT_DIR, "gamma_optical_frequency_vs_nconfigs.dat"),
    output_data,
    header="N_CONFIGS   1/sqrt(N_CONFIGS)   lowest_non_acoustic_Gamma_frequency",
)

print("\nDone.")
print(f"Results saved in: {OUTPUT_DIR}")
