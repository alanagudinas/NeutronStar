############################################################################

#   Minimze difference in mass and radius values for two EOS's to find best fit
#   Author: Alana Gudinas

############################################################################

# Script to find gamma coefficients for continuous equation of state.
# Uses curve-fitting method to minimize difference in values for mass and radius
# between a test equation of state and the reference coefficients for gamma from
# the 7th-order fit (see PP_Gamma_Fit).

############################################################################

import scipy
import numpy as np
from math import log, exp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import TOVSolver_MR_fit as tov
from PresEpsCompMR import CompPresEps
import TOVSolver_differentEOS as tvg
from PP_Gamma_Fit import *
from NeutronStarEOSlibrary import *

# Compute reference vectors for mass and radius values from 7th-order fit to gamma.
# mass0, radius0 = tvg.MRvector0()

# Define vectors of mass and radius results from Piecewise Polytrope with range of central densities:

def rms_error(true,fit):
    return np.sqrt(((true - fit) ** 2).mean())

def MRVecGen(coefvec):
    mass,radius = tov.main(coefvec)
    return mass,radius

def DiffFunc(X,g2,g3):                     # Function of the difference between computed mass/radius and reference mass/radius.
    gamvec = np.array([1.35692,0.0,g2,g3])
    mass, radius = MRVecGen(gamvec)
    return np.sqrt((mass - X[0])**2 + (radius - X[1])**2)

############################################################################
    
xdata = mass0
zdata = radius0
nd = len(xdata)
ydata = np.zeros(nd)
Xdata = zip(xdata,zdata)

# Function that uses curve_fit to minimize DiffFunc. 

def NewMRfit():
    values,covar = curve_fit(DiffFunc,(xdata, zdata), ydata, p0=(0.13641631,-0.01519344))
    return values, covar


# Initial conditions from 5th order fit to gamma:
# p0=(-0.66116166,0.53195721,-0.11793102,0.00811668)

############################################################################
    
# Run curve fit.
    
NewMRfit()

rms_eR = rms_error(mass0,mass)
rms_eM = rms_error(radius0,radius)

print rms_eR, rms_eM
  

