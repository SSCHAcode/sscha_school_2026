import os
import numpy as np
import matplotlib.pyplot as plt

import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.Methods


# ============================================================
# INPUT PARAMETERS
# ============================================================

NQIRR = 3
PATH = "GXMGR"
N_POINTS = 1000

TEMPERATURES = [0,100, 200, 300]

SPECIAL_POINTS = {
    "G": [0.0, 0.0, 0.0],
    "X": [0.5, 0.0, 0.0],
    "M": [0.5, 0.5, 0.0],
    "R": [0.5, 0.5, 0.5],
}

OUTPUT_PLOT = "hessian_dispersion_vs_temperature.png"


# ============================================================
# LOAD REFERENCE STRUCTURE
# ============================================================

first_temperature = TEMPERATURES[0]
first_dyn_prefix = os.path.join(
    f"HESSIAN/T_{first_temperature}",
    "hessian_dyn_"
)

hessian_dyn_ref = CC.Phonons.Phonons(first_dyn_prefix, NQIRR)
nmodes = hessian_dyn_ref.structure.N_atoms * 3


# ============================================================
# GET BAND PATH
# ============================================================

qpath, data = CC.Methods.get_bandpath(
    hessian_dyn_ref.structure.unit_cell,
    PATH,
    SPECIAL_POINTS,
    N_POINTS
)

xaxis, xticks, xlabels = data


# ============================================================
# COLORS
# ============================================================

cmap = plt.get_cmap("coolwarm")
colors = cmap(np.linspace(0, 1, len(TEMPERATURES)))

# ============================================================
# PLOT HESSIAN DISPERSIONS
# ============================================================

plt.figure(dpi=150)
ax = plt.gca()

for itemp, TEMPERATURE in enumerate(TEMPERATURES):

    dyn_prefix = os.path.join(
        f"HESSIAN/T_{TEMPERATURE}",
        "hessian_dyn_"
    )

    if not os.path.exists(dyn_prefix + "1"):
        raise FileNotFoundError(f"Missing file: {dyn_prefix + '1'}")

    print(f"Reading Hessian dynamical matrix: {dyn_prefix}")

    hessian_dyn = CC.Phonons.Phonons(dyn_prefix, NQIRR)
    hessian_dyn.Symmetrize()

    hessian_dispersion = CC.ForceTensor.get_phonons_in_qpath(
        hessian_dyn,
        qpath
    )

    color = colors[itemp]

    for imode in range(nmodes):

        label = None

        if imode == 0:
            label = f"T = {TEMPERATURE} K"

        ax.plot(
            xaxis,
            hessian_dispersion[:, imode],
            color=color,
            lw=1.5,
            label=label
        )


# ============================================================
# AXES STYLE
# ============================================================

ax.legend()

for x in xticks:
    ax.axvline(x, 0, 1, color="k", lw=0.4)

ax.axhline(0, 0, 1, color="k", ls=":", lw=0.4)

ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)

ax.set_xlabel("Q path")
ax.set_ylabel(r"$\omega$ [cm$^{-1}$]")
ax.set_title("Free-energy Hessian dispersion vs temperature")

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300)
plt.show()

print(f"\nPlot saved in: {OUTPUT_PLOT}")
