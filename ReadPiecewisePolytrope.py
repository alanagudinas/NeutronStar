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
from scipy.optimize import curve_fit

from NeutronStarEOSlibrary import *

# G represents capital gamma as specified in Read paper.
# K represents capital kappa.

for i in range(0,3):
    epsilon[i] = ((3e10)**2 + a[i])*rho[i] + K[i]*rho[i]**G[i]/(G[i]-1) 
    a[i+1] = (epsilon[i]/rho[i]) - (3e10)**2 - (K[i+1]/(G[i+1]-1))*rho[i]**(G[i+1]-1)/(3e10)**2
    print a

def p(rhop):

	if 0 < rhop < rho[0]:
		return K[0]*rhop**G[0]

	elif rho[0] < rhop <= rho1:
		return K[1]*rhop**G[1]

	elif rho1 < rhop < rho2:
		return K[2]*rhop**G[2]

	else:
		return K[3]*rhop**G[3]

def epsilon(rhop):
 	if 0 < rhop < rho[0]:
		return K[0]*rhop**G[0]/(G[0]-1) + rhop

	elif rho[0] < rhop <= rho1:
		return  ((3e10)**2+a[1])*rhop + K[1]*rhop**G[1]/(G[1]-1)

	elif rho1 < rhop < rho2:
		return  ((3e10)**2+a[2])*rhop + K[2]*rhop**G[2]/(G[2]-1)

	else:
		return  ((3e10)**2+a[3])*rhop + K[3]*rhop**G[3]/(G[3]-1)



def Gam(rhop):
    if 0 < rhop < rho[0]:
        
        return G[0]
    elif rho[0] < rhop <= rho1:
        
        return G[1]
    elif rho1 < rhop < rho2:
        
        return G[2]
    else:
        return G[3]
    
    
def main(density):
    return epsilon(density)
    return p(density)

    
# Generate vectors of computed values from above functions.
    
LogRhoVec = np.arange(13.,15.,0.01)

for i in range(0,Npts):
    Rho[i] = 10**(LogRhoVec[i])
    Pres[i] = p(Rho[i])
    Epsilon[i] = epsilon(Rho[i])
    Pres[i] = log(Pres[i])
    Rho_log[i] = log(Rho[i])
    Gam_data[i] = Gam(Rho[i])
    xdata[i] = log(Rho[i]/6.1e12)


