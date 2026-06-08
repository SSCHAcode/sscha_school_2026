import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# ============================================================
# INPUT PARAMETERS
# ============================================================

DATA_FILE = "freq_to_plot.dat"

PLOT_FREQ = "hessian_frequency_vs_temperature_fit.png"
PLOT_W2 = "hessian_signed_w2_vs_temperature_fit.png"


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

def signed_sqrt_fit(T, A, Tc):
    return np.sign(T - Tc) * A * np.sqrt(np.abs(T - Tc))


def linear_fit(T, m, q):
    return m * T + q


# ============================================================
# LOAD DATA
# ============================================================

data = np.loadtxt(DATA_FILE)

T = data[:, 0]
omega = data[:, 1]

signed_w2 = np.sign(omega) * omega**2


# ============================================================
# SELECT FIT WINDOWS
# ============================================================

freq_fit_mask = (T >= FIT_FREQ_TMIN) & (T <= FIT_FREQ_TMAX)
w2_fit_mask = (T >= FIT_W2_TMIN) & (T <= FIT_W2_TMAX)

T_freq_fit = T[freq_fit_mask]
omega_fit = omega[freq_fit_mask]

T_w2_fit = T[w2_fit_mask]
signed_w2_fit = signed_w2[w2_fit_mask]

if len(T_freq_fit) < 3:
    raise ValueError("Not enough points in frequency fit window.")

if len(T_w2_fit) < 2:
    raise ValueError("Not enough points in signed-w2 fit window.")


# ============================================================
# FIT 1: signed frequency
# omega(T) = sign(T - Tc) A sqrt(|T - Tc|)
# ============================================================

popt_sqrt, pcov_sqrt = curve_fit(
    signed_sqrt_fit,
    T_freq_fit,
    omega_fit,
    p0=[2.0, 160.0],
    bounds=([0.0, -1000.0], [np.inf, 1000.0])
)

A_fit, Tc_fit = popt_sqrt

T_smooth = np.linspace(np.min(T), np.max(T), 500)
omega_sqrt_fit = signed_sqrt_fit(T_smooth, A_fit, Tc_fit)


# ============================================================
# FIT 2: signed omega squared
# sign(omega) omega^2 = m T + q
# ============================================================

popt_linear, pcov_linear = curve_fit(
    linear_fit,
    T_w2_fit,
    signed_w2_fit
)

m_fit, q_fit = popt_linear
Tc_linear = -q_fit / m_fit

signed_w2_linear_fit = linear_fit(T_smooth, m_fit, q_fit)


# ============================================================
# PLOT 1: omega vs temperature
# ============================================================

plt.figure()

plt.plot(T, omega, "o-", label=r"$\omega$")

plt.plot(
    T_smooth,
    omega_sqrt_fit,
    "--",
    color="green",
    label=(
        rf"Fit $[{FIT_FREQ_TMIN}, {FIT_FREQ_TMAX}]$ K: "
        rf"$T_c={Tc_fit:.2f}$ K"
    )
)

plt.axvspan(FIT_FREQ_TMIN, FIT_FREQ_TMAX, alpha=0.15, label="fit window")
plt.axhline(0, linestyle="--", linewidth=1)

plt.xlabel("Temperature (K)")
plt.ylabel(r"$\omega$ (cm$^{-1}$)")
plt.title(r"Soft-mode frequency $\omega$ vs temperature")

plt.legend()
plt.tight_layout()
plt.savefig(PLOT_FREQ, dpi=300)


# ============================================================
# PLOT 2: sign(omega) omega^2 vs temperature
# ============================================================

plt.figure()

plt.plot(T, signed_w2, "o-", label=r"Eigenvalue")

plt.plot(
    T_smooth,
    signed_w2_linear_fit,
    "--",
    color="green",
    label=(
        rf"Linear fit $[{FIT_W2_TMIN}, {FIT_W2_TMAX}]$ K: "
        rf"$T_c={Tc_linear:.2f}$ K"
    )
)

plt.axvspan(FIT_W2_TMIN, FIT_W2_TMAX, alpha=0.15, label="fit window")
plt.axhline(0, linestyle="--", linewidth=1)

plt.xlabel("Temperature (K)")
plt.ylabel(r"Eigenvalue (cm$^{-1}$)$^2$")
plt.title(r"Soft eigenvalue vs temperature")

plt.legend()
plt.tight_layout()
plt.savefig(PLOT_W2, dpi=300)

plt.show()


# ============================================================
# PRINT FIT RESULTS
# ============================================================

print("\nDone.")
print(f"Data loaded from: {DATA_FILE}")
print(f"Plot saved in: {PLOT_FREQ}")
print(f"Plot saved in: {PLOT_W2}")

print("\nFit windows:")
print(f"Frequency fit window: {FIT_FREQ_TMIN} K <= T <= {FIT_FREQ_TMAX} K")
print(f"Signed-w2 fit window: {FIT_W2_TMIN} K <= T <= {FIT_W2_TMAX} K")

print("\nSigned square-root fit:")
print("omega(T) = sign(T - Tc) A sqrt(|T - Tc|)")
print(f"A  = {A_fit:.8f} cm^-1 K^-1/2")
print(f"Tc = {Tc_fit:.8f} K")

print("\nLinear signed-w2 fit:")
print("sign(omega) omega^2 = m T + q")
print(f"m  = {m_fit:.8f} (cm^-1)^2 / K")
print(f"q  = {q_fit:.8f} (cm^-1)^2")
print(f"Tc = {Tc_linear:.8f} K")
