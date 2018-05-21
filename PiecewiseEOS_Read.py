############################################################################

#   Piecewise Polytrope equation of state from Read paper
#   Author: Alana Gudinas

############################################################################

# Defines functions to compute equation of state as per Read paper.

############################################################################

import scipy
import numpy as np
from math import log, exp
import matplotlib.pyplot as plt
from NeutronStarEOSlibrary import *

## Code unit conversions
M_g = 2e33
M_kg = 2e30
G_c = 6.67e-8
c_c = 3e10
L_cc = (G_c*M_g)/c_c**2
T_cc = (M_g*G_c)/c_c**3

a = np.zeros(4)
epsilon = np.zeros(3) 
p = np.zeros(3)
K = np.zeros(4)

#SLy values for gamma
G = np.array([1.35692, 3.005, 2.988, 2.851])

K[0] = 0.0894844158847432
p[1] = (10**34.4)/((M_g)/((T_cc**2)*L_cc))
rho = np.array([0,8.1604e-4, 1.6282e-3])
rho1 = 8.1604e-4
rho2 = 1.6282e-3

## Calculating values for piecewise polytrope

K[1] = p[1]/(rho1**G[1])
rho[0] = (K[0]/K[1])**(1/(G[1]-G[0]))
p[0] = K[0] * rho[0]**G[0]
K[2] = p[1]/(rho1**G[2])
p[2] = K[2] * rho2**G[2]
K[3] = p[2]/(rho2**G[3])


for i in range(0,3):
    epsilon[i] = ((1)**2 + a[i])*rho[i] + K[i]*rho[i]**G[i]/(G[i]-1) 
    a[i+1] = (epsilon[i]/rho[i]) - (1)**2 - (K[i+1]/(G[i+1]-1))*rho[i]**(G[i+1]-1)/(1)**2
    #print a

def PressureFromDensityGammaPP(rhop):
    if 0 < rhop < rho[0]:
        return K[0]*rhop**G[0]
    elif rho[0] < rhop <= rho1:
        return K[1]*rhop**G[1]

    elif rho1 < rhop < rho2:
        return K[2]*rhop**G[2]

    else:
        return K[3]*rhop**G[3]
    
def EpsilonFromDensityGammaPP(rhop):
    if 0 < rhop < rho[0]:
        return K[0]*rhop**G[0]/(G[0]-1) + rhop

    elif rho[0] < rhop <= rho1:
        return  ((1)**2+a[1])*rhop + K[1]*rhop**G[1]/(G[1]-1)

    elif rho1 < rhop < rho2:
        return  ((1)**2+a[2])*rhop + K[2]*rhop**G[2]/(G[2]-1)

    else:
        return  ((1)**2+a[3])*rhop + K[3]*rhop**G[3]/(G[3]-1)
    
def Gam(rhop):
    if 0 < rhop < rho[0]:
        return G[0]
    elif rho[0] < rhop <= rho1:
        
        return G[1]
    elif rho1 < rhop < rho2:
        
        return G[2]
    else:
        return G[3]
