import cellconstructor as CC, cellconstructor.Phonons
import sys, os

DYN_NAME = "sscha_converged_dyn_"
NQIRR = 4

def extract_structure():
    dyn = CC.Phonons.Phonons(DYN_NAME, NQIRR)
    dyn.structure.save_scf("sscha_structure.scf")


if __name__ == "__main__":
    # Change the working directory to the one containing this script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    extract_structure()
    
