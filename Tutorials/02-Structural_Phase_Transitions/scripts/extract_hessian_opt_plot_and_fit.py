import os
import numpy as np
import matplotlib.pyplot as plt

import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.Units as Units

from scipy.optimize import curve_fit


# ============================================================
# INPUT PARAMETERS
# ============================================================

TEMPERATURES = [0, 50, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300]

DATA_FILE = "freq_to_plot.dat"

PLOT_FREQ = "freq_vs_temperature.png"
PLOT_HESSIAN_W2 = "hessian_eigenvalue_vs_temperature.png"


# ============================================================
# FIT WINDOWS
# ============================================================

FIT_FREQ_TMIN = 100
FIT_FREQ_TMAX = 200

FIT_W2_TMIN = 100
FIT_W2_TMAX = 200


# ============================================================
# FUNCTIONS
# ============================================================

def temperature_label(T):
    return str(T).replace(".", "p")


def get_soft_optical_frequency_cm(dyn_file):
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

    return optical_freqs[np.argmin(np.abs(optical_freqs))]


def signed_sqrt_fit(T, A, Tc):
    return np.sign(T - Tc) * A * np.sqrt(np.abs(T - Tc))


def linear_fit(T, m, q):
    return m * T + q


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
# PREPARE DATA
# ============================================================

T = results[:, 0]
hessian_freq = results[:, 1]
sscha_freq = results[:, 2]

hessian_signed_w2 = np.sign(hessian_freq) * hessian_freq**2


# ============================================================
# SELECT FIT WINDOWS
# ============================================================

freq_fit_mask = (
    (T >= FIT_FREQ_TMIN)
    & (T <= FIT_FREQ_TMAX)
)

w2_fit_mask = (
    (T >= FIT_W2_TMIN)
    & (T <= FIT_W2_TMAX)
)

T_freq_fit = T[freq_fit_mask]
hessian_freq_fit = hessian_freq[freq_fit_mask]

T_w2_fit = T[w2_fit_mask]
hessian_w2_fit = hessian_signed_w2[w2_fit_mask]

if len(T_freq_fit) < 3:
    raise ValueError(
        "Not enough points in frequency fit window. "
        "Increase FIT_FREQ_TMIN/FIT_FREQ_TMAX range."
    )

if len(T_w2_fit) < 2:
    raise ValueError(
        "Not enough points in eigenvalue fit window. "
        "Increase FIT_W2_TMIN/FIT_W2_TMAX range."
    )


# ============================================================
# FIT 1: signed Hessian frequency
# omega(T) = sign(T - Tc) A sqrt(|T - Tc|)
# ============================================================

popt_sqrt, pcov_sqrt = curve_fit(
    signed_sqrt_fit,
    T_freq_fit,
    hessian_freq_fit,
    p0=[2.0, 160.0],
    bounds=([0.0, -1000.0], [np.inf, 1000.0])
)

A_fit, Tc_fit = popt_sqrt

T_smooth_freq = np.linspace(np.min(T), np.max(T), 500)
hessian_freq_sqrt_fit = signed_sqrt_fit(
    T_smooth_freq,
    A_fit,
    Tc_fit
)


# ============================================================
# FIT 2: signed Hessian eigenvalue
# omega^2_signed(T) = m T + q
# ============================================================

popt_linear, pcov_linear = curve_fit(
    linear_fit,
    T_w2_fit,
    hessian_w2_fit
)

m_fit, q_fit = popt_linear

Tc_linear = -q_fit / m_fit

T_smooth_w2 = np.linspace(np.min(T), np.max(T), 500)
hessian_w2_linear_fit = linear_fit(
    T_smooth_w2,
    m_fit,
    q_fit
)


# ============================================================
# PLOT 1: Hessian and SSCHA frequencies vs temperature
# ============================================================

plt.figure()

plt.plot(T, hessian_freq, "o-", label="Hessian frequency")
#plt.plot(T, sscha_freq, "s-", label="SSCHA")

plt.plot(
    T_smooth_freq,
    hessian_freq_sqrt_fit,
    "--",
    color="green",
    label=(
        rf"Hessian fit "
        rf"$[{FIT_FREQ_TMIN}, {FIT_FREQ_TMAX}]$ K: "
        rf"$T_c={Tc_fit:.2f}$ K"
    )
)

plt.axvspan(
    FIT_FREQ_TMIN,
    FIT_FREQ_TMAX,
    alpha=0.15,
    label="fit window"
)

plt.axhline(0, linestyle="--", linewidth=1)

plt.xlabel("Temperature (K)")
plt.ylabel("Optical frequency (cm$^{-1}$)")
plt.title("Optical mode frequency vs temperature")

plt.legend()
plt.tight_layout()
plt.savefig(PLOT_FREQ, dpi=300)


# ============================================================
# PLOT 2: signed Hessian frequency squared vs temperature
# ============================================================

plt.figure()

plt.plot(T, hessian_signed_w2, "o-", label="Hessian eigenvalue")

plt.plot(
    T_smooth_w2,
    hessian_w2_linear_fit,
    "--",
    color="green",
    label=(
        rf"Linear fit "
        rf"$[{FIT_W2_TMIN}, {FIT_W2_TMAX}]$ K: "
        rf"$T_c={Tc_linear:.2f}$ K"
    )
)

plt.axvspan(
    FIT_W2_TMIN,
    FIT_W2_TMAX,
    alpha=0.15,
    label="fit window"
)

plt.axhline(0, linestyle="--", linewidth=1)

plt.xlabel("Temperature (K)")
plt.ylabel(r"Eigenvalue (cm$^{-1}$)$^2$")
plt.title('Hessian "optical eigenvalue" vs temperature')

plt.legend()
plt.tight_layout()
plt.savefig(PLOT_HESSIAN_W2, dpi=300)

plt.show()


# ============================================================
# PRINT FIT RESULTS
# ============================================================

print("\nDone.")
print(f"Data saved in: {DATA_FILE}")
print(f"Plot saved in: {PLOT_FREQ}")
print(f"Plot saved in: {PLOT_HESSIAN_W2}")

print("\nFit windows:")
print(f"Frequency fit window:  {FIT_FREQ_TMIN} K  <= T <= {FIT_FREQ_TMAX} K")
print(f"Eigenvalue fit window: {FIT_W2_TMIN} K  <= T <= {FIT_W2_TMAX} K")

print("\nSigned square-root fit:")
print("omega(T) = sign(T - Tc) A sqrt(|T - Tc|)")
print(f"A  = {A_fit:.8f} cm^-1 K^-1/2")
print(f"Tc = {Tc_fit:.8f} K")

print("\nLinear eigenvalue fit:")
print("omega_signed^2(T) = m T + q")
print(f"m  = {m_fit:.8f} (cm^-1)^2 / K")
print(f"q  = {q_fit:.8f} (cm^-1)^2")
print(f"Tc = {Tc_linear:.8f} K")
