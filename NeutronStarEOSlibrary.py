############################################################################

#   Neutron Star EOS Library
#   Author: Alana Gudinas

############################################################################


# This is a libary containing numerical values needed to compute different 
# equations of state for neutron stars. Sections are separated according
# to different equations of states from various papers.

import scipy
import numpy as np
from math import log, exp


# Lindblom, "Spectral Representations of Neutron-Star Equations of State," 2010.
# Data used in computing values for energy density and pressure (SLy).

p0 = 1.64e33   # Initial pressure in cgs units
g0 = 0.9865        # Spectral coefficients 
g1 = 0.1110
g2 = -0.0301
g3 = 0.0022
eps0 = 2.05e14*(3e10)**2 # Intial energy density

############################################################################

# Read, "Constraints on a phenomenologically parameterized neutron-star equation of state," 2008.
# Data used in computing piecewise polytropes with four parameters. 

c = 1 

# Empty vectors for polytrope variables.

a = np.zeros(4)
epsilon = np.zeros(3)
p = np.zeros(3)
K = np.zeros(4)

# Adiabitic indicies for each polytrope (SLy EOS):

G = np.array([1.3569, 3.005, 2.988, 2.851])

# Coefficient for first polytrope.

K[0] = (3.998736e-8)*(3e10)**2

# Value for p1 (pressure):

p[1] = 10**(34.384)/(c**2)

# Using values for rho1 and rho2 (density) specified in Read paper, compute values for each polytrope:

rho = np.array([0, 10**14.7, 10**15]) # cgs units.

rho1 = 10**14.7
rho2 = 10**15

epsilon[0] = rho[0] + p[0]/(G[0]-1) # first polytrope.

# Equation for pressure as function of density as specified in paper: p(rho) = K*rho^G

K[1] = p[1]/(rho1**G[1])
rho[0] = (K[0]/K[1])**(1/(G[1]-G[0]))
p[0] = K[0] * rho[0]**G[0]
K[2] = p[1]/(rho1**G[2])
p[2] = K[2] * rho2**G[2]
K[3] = p[2]/(rho2**G[3])

# Vectors for plotting purposes:

LogRhoVec = np.arange(13.,15.,0.01)
Npts = len(LogRhoVec)
Rho = np.zeros(Npts)
Pres = np.zeros(Npts)
Epsilon = np.zeros(Npts)
Gam_data = np.zeros(Npts)
xdata = np.zeros(Npts)
Rho_log = np.zeros(Npts)

############################################################################

# Douchin, "A unified equation of state of dense matter and neutron star structure," 2001.

# Vectors of nuclear data for gamma coefficients and densities from crust to core. 

Gamma_nuclear = np.array([1.344,1.353,1.351,1.342,1.332,1.322,1.320,1.325,1.338,1.358,1.387,1.416,1.458,1.496,1.536,1.576,1.615,1.650,1.672,1.686,1.685,1.662,1.644,2.159,2.217,2.309,2.394,2.539,2.655,2.708,2.746,2.905,2.990,3.025,3.035,3.032,3.023,3.012,2.999,2.987,2.975,2.964,2.953,2.943,2.933,2.924,2.916,2.908,2.900,2.893,2.881,2.869,2.858,2.847,2.836,2.824,2.801,2.778,2.754,2.731,2.708])
density_nuclear = np.array([7.5106e12,9.6148e12,1.2593e13,1.6774e13,2.1042e13,2.7844e13,3.6043e13,4.0688e13,4.7001e13,5.3843e13,6.1153e13,6.7284e13,7.5224e13,8.1738e13,8.8350e13,9.5022e13,1.0173e14,1.0845e14,1.1351e14,1.1859e14,1.2372e14,1.272e14,1.2845e14,1.3038e14,1.3531e14,1.4381e14,1.5232e14,1.6935e14,1.8641e14,2.0350e14,2.2063e14,2.7223e14,3.2424e14,3.7675e14,4.2983e14,4.8358e14,5.3808e14,5.9340e14,6.4963e14,7.0684e14,7.6510e14,8.2450e14,8.8509e14,9.4695e14,1.0102e15,1.0748e15,1.1408e15,1.2085e15,1.2777e15,1.3486e15,1.4706e15,1.5977e15,1.7302e15,1.8683e15,2.0123e15,2.1624e15,2.4820e15,2.8289e15,3.2048e15,3.6113e15,4.0498e15])
pressure_nuclear = np.array([1.0405e31, 1.4513e31, 2.0894e31, 3.0720e31, 4.1574e31, 6.0234e31, 8.4613e31, 9.9286e31, 1.2023e32, 1.4430e32, 1.7175e32, 1.9626e32, 2.3024e32, 2.6018e32, 2.9261e32, 3.2756e32, 3.6505e32, 4.0509e32, 4.3681e32, 4.6998e32, 5.0462e32, 5.2856e32, 5.3739e32, 5.3739e32, 5.8260e32, 6.6828e32, 7.6443e32, 9.9146e32, 1.2701e33, 1.6063e33, 1.9971e33, 3.5927e33, 5.9667e33, 9.2766e33, 1.3668e34, 1.9277e34, 2.6235e34, 3.4670e34, 4.4702e34, 5.6451e34, 7.0033e34, 8.5561e34, 1.0315e35, 1.2289e35, 1.4491e35, 1.6930e35, 1.9616e35, 2.2559e35, 2.5769e35, 2.9255e35, 3.5702e35, 4.2981e35, 5.1129e35, 6.0183e35, 7.0176e35, 8.1139e35, 1.0609e36, 1.3524e36, 1.6876e36, 2.0679e36, 2.4947e36])

nd = len(density_nuclear)

G0 = 1.35692    # G0 that aligns with Read and Lindblom papers.

# Conversion factors for code units.

M_g = 2e33
M_kg = 2e30
G_c = 6.67e-8
c_c = 3e10
L_cc = (G_c*M_g)/c_c**2
T_cc = (M_g*G_c)/c_c**3

## Convert nuclear data from cgs to code units.

PresNucCode = pressure_nuclear/((M_g)/((T_cc**2)*L_cc))
DensNucCode = density_nuclear*(L_cc**3)/M_g

DensNucCode_x = np.zeros(nd)
Gam_dataPP = np.zeros(nd)

for i in range(nd):
    DensNucCode_x[i] = log(DensNucCode[i]/DensNucCode[0])

rhostart = 0.0003

############################################################################

# Vectors of gamma coefficients to test various EOS fits in the TOV Solver.

# 7th order fit to piecewise polytrope.

originalGam = np.array([-0.00344500377,0.0725210426,-.579092128,2.14785242,-3.58193674,2.09899986,0,1.35692])
originalGam = originalGam[::-1]

# 5th order fit to piecewise poltryope.

testco5th = np.array([0.00811668,-0.11793102,0.53195721,-0.66116166,0,1.35692])
testco5th = testco5th[::-1]

testco6th = np.array([1.344,0,-0.16073642,0,0.08375492,-0.02445526,0.00190295])

############################################################################

# Values for computing fit to gamma coefficients of piecewise polytrope.

initialPressure = 1.9307289822e-08
initialEpsilon = 1.22829154377e-05


## code units of P0 per Read paper
# initialPressure = 1.05499114846e-06
# initalEpsilon = 0.000236236240043
