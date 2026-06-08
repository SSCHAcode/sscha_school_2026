import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# LOAD DATA
# ============================================================

data = np.loadtxt(
    "HESSIAN_CONVERGENCE_T_0/gamma_optical_frequency_vs_nconfigs.dat"
)

# Columns
N_configs = data[:, 0]
inv_sqrt_N = data[:, 1]
freq = data[:, 2]


# ============================================================
# PLOT: frequency vs N_configs
# ============================================================

plt.figure()

plt.plot(
    N_configs,
    freq,
    "o-",
)

plt.xlabel("Number of configurations")
plt.ylabel("Lowest optical Γ frequency (cm-1)")
plt.title("Convergence with ensemble size")

plt.tight_layout()

plt.savefig("freq_vs_N.png", dpi=300)


# ============================================================
# PLOT: frequency vs 1/sqrt(N)
# ============================================================

plt.figure()

plt.plot(
    inv_sqrt_N,
    freq,
    "o-",
)

plt.xlabel(r"$1/\sqrt{N_{\mathrm{configs}}}$")
plt.ylabel("Lowest optical Γ frequency (cm-1)")
plt.title("Monte Carlo convergence scaling")

plt.tight_layout()

plt.savefig("freq_vs_inv_sqrtN.png", dpi=300)

plt.show()
