import cellconstructor as CC, cellconstructor.Units
import tdscha, tdscha.DynamicalLanczos
import numpy as np
import matplotlib.pyplot as plt
import sys, os

W_START = 0 # cm-1
W_END = 200 # cm-1
TERMINATOR = True
SMEARING = 0.5 # cm-1
TEMPERATURE = 300 # K

def plot_spectrum():
    # Get the frequency
    w_array = np.linspace(W_START, W_END, 2000)
    w_ry = w_array / CC.Units.RY_TO_CM
    smearing_ry = SMEARING / CC.Units.RY_TO_CM

    unpol_spectrum = np.zeros_like(w_array)
    
    for i in range(7):
        # Load the final lanczos results
        lanczos = tdscha.DynamicalLanczos.Lanczos()
        lanczos.load_status(f"Raman_unpol{i}.npz")

        green_function = lanczos.get_green_function_continued_fraction(w_ry, smearing=smearing_ry,
                                                                       use_terminator = TERMINATOR)

        # The IR spectrum is proportional to - Im (G(w))
        factor = 45
        if i > 0:
            factor = 7

        unpol_spectrum += -np.imag(green_function) * factor

    # Add the Bose-occupation scattering factor (Stokes)
    unpol_spectrum *= tdscha.DynamicalLanczos.bose_occupation(w_ry, TEMPERATURE) + 1.0

    # Plot the spectrum
    plt.plot(w_array, unpol_spectrum)
    plt.xlabel("Frequency [cm-1]")
    plt.ylabel(r"- Im $(G) (n(\omega) + 1)$")
    plt.title("Raman spectrum - unpolarized")
    plt.tight_layout()
    plt.savefig("raman_spectrum.png")
    plt.show()


if __name__ == "__main__":
    # Change the working directory to the one containing this script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    plot_spectrum()
