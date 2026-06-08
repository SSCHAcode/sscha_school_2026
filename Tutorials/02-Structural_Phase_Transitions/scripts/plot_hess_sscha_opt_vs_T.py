import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("data_to_plot.dat", comments="#")

T = data[:, 0]
hessian_freq = data[:, 1]
sscha_freq = data[:, 2]

# Signed squared Hessian frequency
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
plt.savefig("freq_vs_temperature.png", dpi=300)


# ============================================================
# Plot 2: signed Hessian frequency squared vs temperature
# ============================================================

plt.figure()

plt.plot(T, hessian_signed_w2, "o-")
plt.axhline(0, linestyle="--", linewidth=1)

plt.xlabel("Temperature (K)")
plt.ylabel(r"Eigenvalue (cm$^{-1})^2$")
plt.title("Hessian ''optical eigenvalue'' vs temperature")

plt.tight_layout()
plt.savefig("hessian_eigenvalue_vs_temperature.png", dpi=300)

plt.show()
