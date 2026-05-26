from __future__ import print_function
from __future__ import division
import sys,os

import warnings

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

# Import matplotlib and numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib import cm
import numpy as np

# NumPy moved ComplexWarning in newer versions. This block keeps the script
# compatible with both older and newer NumPy releases.
try:
    ComplexWarning = np.exceptions.ComplexWarning
except AttributeError:
    ComplexWarning = np.ComplexWarning

warnings.filterwarnings("ignore", category=ComplexWarning)

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
dyn = CC.Phonons.Phonons(START_DYN, NQIRR)

# We dertermine the potential to be used
POTENTIAL ='../../Materials/GAP_1.xml'
calc = Potential("IP GAP", param_filename=POTENTIAL)

# We load the ensemble and read the results of the configurations
ensemble = sscha.Ensemble.Ensemble(dyn, TEMPERATURE)
ensemble.load(folder_with_ensemble, population = POPULATION, N = N_RANDOM)
ensemble.has_stress = True

# Define the minimization

minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# Define the steps for the centroids and the force constants

minimizer.min_step_dyn = 0.005         # The minimization step on the dynamical matrix
minimizer.min_step_struc = 0.05        # The minimization step on the structure
minimizer.kong_liu_ratio = 0.2         # The parameter that estimates whether the ensemble is still good
minimizer.meaningful_factor = 0.000001 # How much small the gradient should be before I stop?

# Let's start the minimization

minimizer.init()
minimizer.run()

# Save the minimization details
ioinfo = sscha.Utilities.IOInfo()
ioinfo.SetupSaving("minim_{}".format(POPULATION))

minimizer.run(custom_function_post = ioinfo.CFP_SaveAll)
minimizer.finalize()
namefile='dyn_end_population'+str(POPULATION)+'_'
minimizer.dyn.save_qe(namefile)
