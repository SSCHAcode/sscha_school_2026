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


PATH_HESSIAN="./HESSIAN_T_0_N_256"
PATH_SSCHA = "./RELAX"

sscha_dyn = CC.Phonons.Phonons(
    f"{PATH_SSCHA}/T_0/sscha_dyn_",
    NQIRR
)

hessian_dyn = CC.Phonons.Phonons(
    f"{PATH_HESSIAN}/hessian_dyn_",
    NQIRR
)

###########################################################################################################
nmodes = sscha_dyn.structure.N_atoms * 3

# Get band path
qpath, data = CC.Methods.get_bandpath(
    sscha_dyn.structure.unit_cell,
    PATH,
    SPECIAL_POINTS,
    N_POINTS)

xaxis, xticks, xlabels = data

# Fourier interpolation
sscha_dispersion = CC.ForceTensor.get_phonons_in_qpath(
    sscha_dyn,
    qpath)

hessian_dispersion = CC.ForceTensor.get_phonons_in_qpath(
    hessian_dyn,
    qpath)

plt.figure(dpi=150)
ax = plt.gca()

for i in range(nmodes):
    
    lbl_hessian = None
    lbl_sscha = None

    if i == 0:
        lbl_hessian = "Hessian"
        lbl_sscha = f"SSCHA"

    # Harmonic
    ax.plot(
        xaxis,
        hessian_dispersion[:, i],
        color="b",
        ls="-",
        lw=1,
        label=lbl_hessian)

    # SSCHA
    ax.plot(
        xaxis,
        sscha_dispersion[:, i],
        color="r",
        ls="--",
        lw=1,
        label=lbl_sscha)

#ax.legend()
ax.legend(
    loc="lower center",
    bbox_to_anchor=(0.5, 0.02),
    ncol=1
)


for x in xticks:
    ax.axvline(x, 0, 1, color="k", lw=0.4)

ax.axhline(0, 0, 1, color="k", ls=":", lw=0.4)

ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)

ax.set_xlabel("Q path")
ax.set_ylabel(r"$\omega$ [cm$^{-1}$]")

ax.set_title(f"Free-energy Hessian vs SSCHA phonon dispersion (T= 0 K)")

plt.tight_layout()
plt.savefig("sscha_vs_hessian.png", dpi=300)
plt.show()
