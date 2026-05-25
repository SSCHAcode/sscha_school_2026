import cellconstructor as CC, cellconstructor.Units
import tdscha, tdscha.DynamicalLanczos
import numpy as np
import matplotlib.pyplot as plt
import sys, os

DATA = "IR_x.npz"
W_START = 0 # cm-1
W_END = 200 # cm-1
TERMINATOR = True
SMEARING = 0.5 # cm-1

def plot_spectrum():
    # Load the final lanczos results
    lanczos = tdscha.DynamicalLanczos.Lanczos()
    lanczos.load_status(DATA)

    # Extract the Green function
    w_array = np.linspace(W_START, W_END, 2000)
    w_ry = w_array / CC.Units.RY_TO_CM
    smearing_ry = SMEARING / CC.Units.RY_TO_CM
    green_function = lanczos.get_green_function_continued_fraction(w_ry, smearing=smearing_ry,
                                                                   use_terminator = TERMINATOR)

    # The IR spectrum is proportional to - Im (G(w))
    ir_spectrum = -np.imag(green_function)

    # Plot the IR spectrum
    plt.plot(w_array, ir_spectrum)
    plt.xlabel("Frequency [cm-1]")
    plt.ylabel("- Im (G)")
    plt.title("IR spectrum - polarization x")
    plt.tight_layout()
    plt.savefig("ir_spectrum.png")
    plt.show()


if __name__ == "__main__":
    # Change the working directory to the one containing this script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    plot_spectrum()
