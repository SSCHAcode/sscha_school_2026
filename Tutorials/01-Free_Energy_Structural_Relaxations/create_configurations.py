# Import cellconstructor needed things
import cellconstructor as CC
import cellconstructor.Phonons

# Import the modules to run the sscha
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

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

# We tell where the starting dynamical matrices in reciprocal space are
# for the q points commensurate with the supercell we want to use for
# the calculation. These dynamical matrices should be in Quantum Espresso
# format.
START_DYN = "../../Materials/CsPbI3_cubic_harmonic_"
# In case we want to create the configurations from the result of a previous minimimization we ca use this starting dynamcial matrices
#START_DYN = namefile='dyn_end_population'+str(POPULATION-1)+'_'
#
# We indicate how many irreducible q points are there for the q points in the grid.
NQIRR = 4

#----------------------------------------------

# Step 1: load the harmonic dynamical matrices and the crystal structure

dyn = CC.Phonons.Phonons(START_DYN, NQIRR) # Load them and read the structure
dyn.ForcePositiveDefinite()                # Force positive phonons in case there are imaginary phonon frequencies
dyn.Symmetrize()                           # Impose symmetries and the acoustic sum rule

# Step 2: create the ensemble, the configurations for which we have to calculate the
# forces, energies and stresses

ensemble = sscha.Ensemble.Ensemble(dyn, TEMPERATURE)
ensemble.generate(N_RANDOM)

# Step 3: We save the ensemble

namefile='population'+str(POPULATION)+'_ensemble' # Name of the folder
ensemble.save(namefile, POPULATION)

