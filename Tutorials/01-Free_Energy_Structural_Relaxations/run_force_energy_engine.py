#
# THIS IS AN EXAMPLE FOR THE SSCHA 2026 SCHOOL
# FOR THE CsPbI3 SYSTEM WITH THE USE OF A MACHINE
# LEARNING POTENTIAL.
#
# THIS TUTORIAL SHOWS HOW TO PERFORM SIMPLE SSCHA RELAXATION
# TO CALCULATE THERMODYNAMICAL PROPERTIES SUCH AS STRUCTURE
# RELAXATIONS
#

# Import cellconstructor needed things
import cellconstructor as CC
import cellconstructor.Phonons

# Import the modules to run the sscha
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

# Import the module to be able to run the ML potential
from quippy.potential import Potential

# Import ASE for reading data
import ase

#import sys, os

#----------------------------------------------
# We set up the parameters for the calculation

# The temperature for the calculation in K
TEMPERATURE = 450

# We tell the system that we are goin to create the configurations
# of the first population. This is just used for labeling purposes.
POPULATION = 1

# We determine the number of configurations that we will generate
# for this population
N_RANDOM = 50

# We define the dynamical matrices that generated the ensemble and the folder where
# the ensemble is stored
START_DYN = 'dyn_start_population'+str(POPULATION)+'_'
folder_with_ensemble = 'population'+str(POPULATION)+'_ensemble'
NQIRR = 4

# We dertermine the potential to be used
POTENTIAL ='../../Materials/GAP_1.xml'
calc = Potential("IP GAP", param_filename=POTENTIAL)

# We need to determine how many irreducible q points are there in
# the grid commemsurate with the suprecell. In this example
#----------------------------------------------

RyToEv = 13.605693
BohrToAngstrom = 0.529177

energy_file=str(folder_with_ensemble)+'/energies_supercell_population'+str(POPULATION)+'.dat'
with open(energy_file, "w") as f_energy:

    struct = CC.Structure.Structure()
    for i in range(N_RANDOM):
        namefile=str(folder_with_ensemble)+'/scf_population'+str(POPULATION)+'_'+str(i+1)+'.dat'
        struct.read_scf(namefile)
        ase_struct = struct.get_ase_atoms()
        ase_struct.set_calculator(calc)
        energy = ase_struct.get_potential_energy()
        forces = ase_struct.get_forces()
        stress = ase_struct.get_stress()

        # --- Write energy in Ry ---
        f_energy.write(f"{energy/RyToEv:.10f}\n")

        # --- Write forces file in Ry/Bohr ---
        forces_file = str(folder_with_ensemble)+'/forces_population'+str(POPULATION)+'_'+str(i+1)+'.dat'
        with open(forces_file, "w") as f_forces:
            for f in forces:
                f_forces.write(f"{f[0]* (BohrToAngstrom / RyToEv):.10f} {f[1]* (BohrToAngstrom / RyToEv):.10f} {f[2]* (BohrToAngstrom / RyToEv):.10f}\n")

        # --- Write stress file in Ry/Bohr^3 ---
        stress_file = str(folder_with_ensemble)+'/pressure_population'+str(POPULATION)+'_'+str(i+1)+'.dat'
        with open(stress_file, "w") as f_stress:
            # ASE gives 6 components: xx yy zz yz xz xy
            f_stress.write(f"{stress[0]* (BohrToAngstrom**3 / RyToEv)}  {stress[5]* (BohrToAngstrom**3 / RyToEv)} {stress[4]* (BohrToAngstrom**3 / RyToEv)}\n")
            f_stress.write(f"{stress[5]* (BohrToAngstrom**3 / RyToEv)}  {stress[1]* (BohrToAngstrom**3 / RyToEv)} {stress[3]* (BohrToAngstrom**3 / RyToEv)}\n")
            f_stress.write(f"{stress[4]* (BohrToAngstrom**3 / RyToEv)}  {stress[3]* (BohrToAngstrom**3 / RyToEv)} {stress[2]* (BohrToAngstrom**3 / RyToEv)}\n")
