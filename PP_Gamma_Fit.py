############################################################################

#   Multiple fits to Pieciewise Polytrope gamma coefficients
#   Author: Alana Gudinas

############################################################################

# Purpose of script is to try different-order polynomials to fit to Piecewise Polytrope 
# in order to have a smooth continuous function to compute equation of state while retaining
# accuracy in computations.


# Uses simple polynomial model to fit to G(x), using the values for Gamma defined in Read paper.
# xdata = ln(rho/rho0), ydata is gamma coefficients.
# Curve-fitting method is scipy's curve_fit.

############################################################################

import scipy
import numpy as np
from math import log, exp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.integrate as integrate

from PiecewiseEOS_Read import *
from NeutronStarEOSlibrary import * 

# Different functions defined for different-order fits.
# Testing accuracy of different numbers of free parameters. 

def piecewise_fit(rhop,a,b,c,d,e,f):
    return a*rhop**7 + b*rhop**6 +  c*rhop**5 + d*rhop**4 + e*rhop**3 + f*rhop**2 + 0*rhop + 1.35692

def piecewise_fit3(rhop,a,b):
    return a*rhop**3 + b*rhop**2 + 0*rhop + 1.35692
    
def rms_error(true,fit):
    return np.sqrt(((true - fit) ** 2).mean())

for i in range(nd):
    Gam_dataPP[i] = Gam(DensNucCode[i])
    
values,covar = curve_fit(piecewise_fit, DensNucCode_x, Gam_dataPP)
values3,covar3 = curve_fit(piecewise_fit3, DensNucCode_x, Gam_dataPP)

############################################################################

# Plotting original piecewise polytrope G(x) with fit.

plt.plot(DensNucCode_x, piecewise_fit3(DensNucCode_x, *values3),label="Fit to Gamma")
plt.plot(DensNucCode_x, Gam_dataPP,label="Gamma")
plt.xlabel("x = log(rho/rho0)")
plt.ylabel("Gamma coefficients")
plt.legend(loc="upper left")
plt.title("Fit to Piecewise Polytrope Gamma")
plt.show()

print values3

rms_e = rms_error(Gam_dataPP, piecewise_fit(DensNucCode_x, *values))

############################################################################

# Computing polynomial using Horner's method.

def polynomial(x,A):
    n = len(A)
    y = A[n-1]*x
    for i in range (n-2,0,-1):
        y = (y + A[i])*x
    return y + A[0]

############################################################################
    
# Vector of coefficients for fit:
    
gamPPcoef = np.array([ -3.44500377e-03, 7.25210426e-02, -5.79092128e-01, 2.14785242e+00, -3.58193674e+00, 2.09899986e+00, 0, 1.35692])

gamPPcoef = gamPPcoef[::-1]

############################################################################

# Compute values of gamma based on different values of x = ln(rho/rho0).

def Gamma_function(x):
    if x < 0:
        return G0
    return polynomial(x,gamPPcoef)

# Integrate Gamma_function to compute pressure and energy density:

def PressureFromGammaPP(rho):
    if rho < 1.e-20:
        return 1.e-20
    x = log(rho/DensNucCode[0])
    Gamma_function_integral= integrate.quad(lambda q: Gamma_function(q), 0, x)[0]
    return initialPressure*exp(Gamma_function_integral)
    
    
def EpsilonFromGammaPP(rho):
    if rho < 1.e-20:
        return 1.e-20
    x = log(rho/DensNucCode[0])
    Pressure_function_integral = integrate.quad(lambda q: PressureFromGammaPP(DensNucCode[0]*exp(q))*exp(-q), 0, x)[0]
    return initialEpsilon*exp(x) + exp(x)*Pressure_function_integral

############################################################################



