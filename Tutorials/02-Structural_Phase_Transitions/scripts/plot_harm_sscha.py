import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.Methods

import numpy as np
import matplotlib.pyplot as plt

NQIRR = 3
PATH = "GXMGR"
N_POINTS = 1000

SPECIAL_POINTS = {
    "G": [0.0, 0.0, 0.0],
    "X": [0.5, 0.0, 0.0],
    "M": [0.5, 0.5, 0.0],
    "R": [0.5, 0.5, 0.5],
}


PATH_HARM="../../Materials/tutorial_02/"
PATH_SSCHA = "./RELAX"

TEMPERATURE = 0

sscha_dyn = CC.Phonons.Phonons(
    f"{PATH_SSCHA}/T_{TEMPERATURE}/sscha_dyn_",
    NQIRR
)

###########################################################################################################
# Load HARM phonons
harm_dyn = CC.Phonons.Phonons(f"{PATH_HARM}ffield_dynq_", NQIRR)
nmodes = harm_dyn.structure.N_atoms * 3

# Get band path
qpath, data = CC.Methods.get_bandpath(
    harm_dyn.structure.unit_cell,
    PATH,
    SPECIAL_POINTS,
    N_POINTS)

xaxis, xticks, xlabels = data

# Fourier interpolation
sscha_dispersion = CC.ForceTensor.get_phonons_in_qpath(
    sscha_dyn,
    qpath)

harm_dispersion = CC.ForceTensor.get_phonons_in_qpath(
    harm_dyn,
    qpath)

plt.figure(dpi=150)
ax = plt.gca()

for i in range(nmodes):
    
    lbl_harm = None
    lbl_sscha = None

    if i == 0:
        lbl_harm = "Harmonic"
        lbl_sscha = f"SSCHA (T={TEMPERATURE} K)"

    # Harmonic
    ax.plot(
        xaxis,
        harm_dispersion[:, i],
        color="k",
        ls="--",
        lw=1,
        label=lbl_harm)

    # SSCHA
    ax.plot(
        xaxis,
        sscha_dispersion[:, i],
        color="r",
        lw=1,
        label=lbl_sscha)

ax.legend()

for x in xticks:
    ax.axvline(x, 0, 1, color="k", lw=0.4)

ax.axhline(0, 0, 1, color="k", ls=":", lw=0.4)

ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)

ax.set_xlabel("Q path")
ax.set_ylabel(r"$\omega$ [cm$^{-1}$]")

plt.tight_layout()
plt.savefig("harm_vs_sscha.png", dpi=300)
plt.show()
