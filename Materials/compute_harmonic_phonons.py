# This simple script computes the harmonic phonons of the given structure.
# We need to use the quippy GAP potential
import sys, os
import cellconstructor as CC, cellconstructor.Phonons

# Load the quippy interatomic potential (ASE calculator)
from quippy.potential import Potential

# Load the interatomic Potential for CsPbI3
current_dir = os.path.dirname(os.path.abspath(__file__))
potential_path = os.path.join(current_dir, "GAP_1.xml")
calc = Potential("IP GAP", param_filename=potential_path)

# Load the CsPbI3 structure
structure_path = os.path.join(current_dir, "CsPbI3_cubic_phase.scf")
structure = CC.Structure.Structure()
structure.read_scf(structure_path)

# We need the atomic masses to compute phonons, which are not in the given file
# ASE provides the masses for every element, so we can use that
structure.build_masses()

# Compute the harmonic phonons using finite displacements
phonons = CC.Phonons.compute_phonons_finite_displacements(structure, calc,
                                                          supercell = (2, 2, 2)
                                                          )

# Save the phonons in quantum-espresso format
phonons_path = os.path.join(current_dir, "CsPbI3_cubic_harmonic_")
phonons.save_qe(phonons_path)
print("Done.")
