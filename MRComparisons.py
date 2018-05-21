############################################################################

#   Comparison of Mass-Radius Curves for different EOS's
#   Author: Alana Gudinas

############################################################################

# Plots various Mass-Radius results.

############################################################################

import scipy
import numpy as np
import matplotlib.pyplot as plt
import TOVSolver_MR_fit as tov
import TOVSolver_differentEOS as tvg
from NeutronStarEOSlibrary import *

# Function that runs TOV Solver for gamma coefficient vector of any length.


def MRComp(coefvec):
    mass,radius = tov.main(coefvec)
    return mass,radius


# Call function for different coefficients to compute mass and radius vectors.

# Original Piecewise Polytrope from Read:

mass0,radius0 = tvg.MRvector0()

# Three different fits to the Piecewise Polytrope, using coefficients already computed:

mass3,radius3 = MRComp(np.array([1.35692, 0, 0.13641631, -0.01519344]))
mass5,radius5 = MRComp(np.array([1.35692, 0, -0.82507687, 0.63569048, -0.13953949, 0.0095992]))
mass7,radius7 = MRComp(np.array([1.35692, 0, 2.09899986, -3.58193674, 2.14785242, -.579092128, .0725210426, -.00344500377]))

# Four different fits to nuclear data:

mass_n3, radius_n3 = MRComp(np.array([G0, -0.65844171,  0.46372075, -0.05281594]))
mass_n5, radius_n5 = MRComp(np.array([G0, 0, -0.65166217, 0.52100071, -0.11468986, 0.00783277 ]))
mass_n7, radius_n7 = MRComp(np.array([G0, 0, 1.10921894, -2.0444276, 1.28847678, -0.35593424, 0.04509961, -0.00215455]))
mass_n15, radius_n15 = MRComp(np.array([G0, 0, -1.36830300,   1.25742819e+01,  -4.04451816e+01,   6.60371042e+01, -6.43776195e+01, 4.09064726e+01, -1.78578819e+01, 5.53075281e+00, -1.23529520e+00, 1.99150632e-01, -2.27538643e-02, 1.75629122e-03, -8.23540471e-05, 1.77279501e-06]))

# Plot all the mass-radius curves together.

plt.plot(radius0,mass0,label="Piecewise Polytrope")
plt.plot(radius3,mass3,label="3rd (PP EOS)")
plt.plot(radius5,mass5,label="5th (PP EOS)")
plt.plot(radius7,mass7,label="7th (PP EOS)")
plt.plot(radius_n3,mass_n3,label="3rd (Nuc EOS)")
plt.plot(radius_n5,mass_n5,label="5th (Nuc EOS)")
plt.plot(radius_n7,mass_n7,label="7th (Nuc EOS)")
plt.plot(radius_n15,mass_n15,label="15th (Nuc EOS)")
plt.xlabel("Radius")
plt.ylabel("Solar Mass")
plt.legend(loc="upper right")
plt.title("Comparison of Mass/Radius for Various EOS's")
plt.grid(True)
plt.show()