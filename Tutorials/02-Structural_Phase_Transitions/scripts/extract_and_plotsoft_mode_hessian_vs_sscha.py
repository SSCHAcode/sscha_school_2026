import os
import numpy as np
import matplotlib.pyplot as plt

import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.Units as Units


# ============================================================
# INPUT PARAMETERS
# ============================================================

TEMPERATURES = [0, 50, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300]

DATA_FILE = "freq_to_plot.dat"

PLOT_FREQ = "freq_vs_temperature.png"
PLOT_HESSIAN_W2 = "hessian_eigenvalue_vs_temperature.png"


# ============================================================
# FUNCTIONS
# ============================================================

def temperature_label(T):
    return str(T).replace(".", "p")


def get_soft_optical_frequency_cm(dyn_file):
    """
    Extract the soft optical frequency in cm^-1 from a single QE dyn file.
    The 3 acoustic modes are removed as the three modes closest to zero.
    """

    if not dyn_file.endswith("_1"):
        raise ValueError(f"Expected file ending with '_1', got: {dyn_file}")

    dyn_prefix = dyn_file[:-1]

    dyn = CC.Phonons.Phonons(dyn_prefix, 1)
    dyn.Symmetrize()

    freqs, pols = dyn.DiagonalizeSupercell()

    freqs_cm = np.real(freqs) * Units.RY_TO_CM

    acoustic_indices = np.argsort(np.abs(freqs_cm))[:3]

    mask = np.ones(len(freqs_cm), dtype=bool)
    mask[acoustic_indices] = False

    optical_freqs = freqs_cm[mask]

    soft_optical_freq = optical_freqs[np.argmin(np.abs(optical_freqs))]

    return soft_optical_freq


# ============================================================
# EXTRACT FREQUENCIES
# ============================================================

results = []

for T in TEMPERATURES:

    tlabel = temperature_label(T)
    folder = f"HESSIAN/T_{tlabel}"

    hessian_file = os.path.join(folder, "hessian_dyn_1")
    sscha_file = os.path.join(folder, "refined_sscha_dyn_1")

    if not os.path.exists(hessian_file):
        raise FileNotFoundError(f"Missing file: {hessian_file}")

    if not os.path.exists(sscha_file):
        raise FileNotFoundError(f"Missing file: {sscha_file}")

    hessian_freq = get_soft_optical_frequency_cm(hessian_file)
    sscha_freq = get_soft_optical_frequency_cm(sscha_file)

    results.append([T, hessian_freq, sscha_freq])

    print(
        f"T = {T:8.3f} K   "
        f"Hessian = {hessian_freq:12.6f} cm^-1   "
        f"SSCHA = {sscha_freq:12.6f} cm^-1"
    )


results = np.array(results)


# ============================================================
# SAVE TABLE
# ============================================================

np.savetxt(
    DATA_FILE,
    results,
    header="Temperature(K)    Hessian(cm^-1)    SSCHA(cm^-1)",
    fmt="%12.6f  %16.8f  %16.8f"
)


# ============================================================
# PLOTS
# ============================================================

T = results[:, 0]
hessian_freq = results[:, 1]
sscha_freq = results[:, 2]

hessian_signed_w2 = np.sign(hessian_freq) * hessian_freq**2


# ============================================================
# Plot 1: Hessian and SSCHA frequencies vs temperature
# ============================================================

plt.figure()

plt.plot(T, hessian_freq, "o-", label="Hessian")
plt.plot(T, sscha_freq, "s-", label="SSCHA")

plt.axhline(0, linestyle="--", linewidth=1)

plt.xlabel("Temperature (K)")
plt.ylabel("Optical frequency (cm$^{-1}$)")
plt.title("Optical mode frequency vs temperature")

plt.legend()
plt.tight_layout()
plt.savefig(PLOT_FREQ, dpi=300)


# ============================================================
# Plot 2: signed Hessian frequency squared vs temperature
# ============================================================

plt.figure()

plt.plot(T, hessian_signed_w2, "o-")
plt.axhline(0, linestyle="--", linewidth=1)

plt.xlabel("Temperature (K)")
plt.ylabel(r"Eigenvalue (cm$^{-1}$)$^2$")
plt.title('Hessian "optical eigenvalue" vs temperature')

plt.tight_layout()
plt.savefig(PLOT_HESSIAN_W2, dpi=300)

plt.show()


print("\nDone.")
print(f"Data saved in: {DATA_FILE}")
print(f"Plot saved in: {PLOT_FREQ}")
print(f"Plot saved in: {PLOT_HESSIAN_W2}")
