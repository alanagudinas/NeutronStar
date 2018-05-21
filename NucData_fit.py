############################################################################

#   Multiple fits to gamma coefficients from nuclear data
#   Author: Alana Gudinas

############################################################################

# Purpose of script is to try different-order polynomials to fit to nuclear data from Douchin paper.
# Will hopefully yield more accurate results for mass and radius computations.

# Uses simple polynomial model to fit to G(x), using the values for Gamma defined in Douchin paper.
# xdata = ln(rho/rho0), ydata is gamma coefficients.
# Curve-fitting method is scipy's curve_fit.

############################################################################

import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from NeutronStarEOSlibrary import *


def rms_error(true,fit):
    return np.sqrt(((true - fit) ** 2).mean())

# Starting at rho = 7.5e12 (cgs units)

# The 15th order polynomial is to achieve a high-accuracy fit to the nuclear data to test in TOV solver. 
    
def nuclear_fit15(rhop,a,b,c,d,e,f,g,h,i,j,k,l,m,n):
    return a*rhop**15 + b*rhop**14 + c*rhop**13 + d*rhop**12 + e*rhop**11 + f*rhop**10 + g*rhop**9 + h*rhop**8 + i*rhop**7 + j*rhop**6 + k*rhop**5 + l*rhop**4 + m*rhop**3 + n*rhop**2 + 0*rhop + G0

def nuclear_fit7(rhop,a,b,c,d,e,f):
    return a*rhop**7 + b*rhop**6 + c*rhop**5 + d*rhop**4 + e*rhop**3 + f*rhop**2 + 0*rhop + G0


def nuclear_fit5(rhop,a,b,c,d):
    return a*rhop**5 + b*rhop**4 + c*rhop**3 + d*rhop**2 + 0*rhop + G0


def nuclear_fit3(rhop,a,b,c):
    return a*rhop**3 + b*rhop**2 + c*rhop + G0

val7,cov7 = curve_fit(nuclear_fit7, DensNucCode_x, Gamma_nuclear)
print val7[::-1]
val5,cov5 = curve_fit(nuclear_fit5, DensNucCode_x, Gamma_nuclear)
print val5[::-1]
val3,cov3 = curve_fit(nuclear_fit3, DensNucCode_x, Gamma_nuclear)
print val3[::-1]
val15,cov15 = curve_fit(nuclear_fit15, DensNucCode_x, Gamma_nuclear)
print val15[::-1]
# Plot original nuclear data along with the four fits:

plt.plot(DensNucCode_x, nuclear_fit7(DensNucCode_x, *val7),label="7th order fit")
plt.plot(DensNucCode_x, nuclear_fit5(DensNucCode_x, *val5),label="5th order fit")
plt.plot(DensNucCode_x, nuclear_fit3(DensNucCode_x, *val3),label="3rd order fit")
plt.plot(DensNucCode_x, nuclear_fit15(DensNucCode_x, *val15),label="15th order fit")
plt.plot(DensNucCode_x, Gamma_nuclear,label="Nuclear Data")
plt.xlabel("x = log(rho/rho0)")
plt.ylabel("Gamma coefficients")
plt.legend(loc="upper left")
plt.title("Fit to Gamma from Nuclear Data")
plt.show()
rms_e = rms_error(Gamma_nuclear, nuclear_fit7(DensNucCode_x, *val7))

print rms_e

