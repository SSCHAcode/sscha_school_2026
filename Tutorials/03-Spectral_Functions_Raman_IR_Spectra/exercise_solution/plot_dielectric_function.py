import cellconstructor as CC, cellconstructor.Units
import tdscha
import tdscha.DynamicalLanczos
import numpy as np
import matplotlib.pyplot as plt

# Load the dynamical matrix
dyn = CC.Phonons.Phonons("sscha_converged_dyn_", 4)

# Load the dielectric tensor
dyn.ReadInfoFromESPRESSO("dielectric_calc.pho")

# Load the tdscha result
lanczos = tdscha.DynamicalLanczos.Lanczos()
lanczos.load_status("IR_x.npz")

# Get the frequencies
W_START = 0
W_END = 300 # cm-1
N_W = 1000
SMEARING = 0.5 # cm-1
w = np.linspace(W_START, W_END, N_W)

# Convert in Ry
w_ry = w / CC.Units.RY_TO_CM
smearing = SMEARING / CC.Units.RY_TO_CM

# Get the suscieptibility
chi = lanczos.get_green_function_continued_fraction(w_ry, 
                                                    smearing=smearing)

# Get the ionic dielectric tensor
volume = dyn.structure.get_volume() * np.prod(dyn.GetSupercell())
volume /= CC.Units.BOHR_TO_ANGSTROM**3 # convert in Borh^3
epsilon_ion = (8*np.pi / volume) * chi

# Get the total dielectric tensor
# Sum the electronic (xx component) and the ionic one
epsilon = dyn.dielectric_tensor[0,0] + epsilon_ion

# Compute the normal propagation and evaneshent wave
# And the standard spectral function
normal_prop = - np.imag(np.sqrt(epsilon))
evaneshent_prop = np.imag(1 / epsilon)
spectral_func = -np.imag(chi)

# Plot the results
fig, axarr = plt.subplots(ncols=1, nrows=3, sharex=True)
axarr[0].plot(w, spectral_func)
axarr[1].plot(w, normal_prop)
axarr[2].plot(w, evaneshent_prop)
axarr[2].set_xlabel("frequency [cm-1]")
axarr[0].set_ylabel(r"- Im $\chi$")
axarr[1].set_ylabel(r"- Im $n$")
axarr[2].set_ylabel(r"- Im $\frac{1}{\epsilon}$")
axarr[0].set_title("Spectral function")
axarr[1].set_title("Normal wave absorbance")
axarr[2].set_title("Evaneshent wave absorbance")
fig.tight_layout()
plt.show()
fig.savefig("dielectric_function.png")
