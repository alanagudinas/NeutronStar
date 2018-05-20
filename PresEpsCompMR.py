############################################################################

#   Compute pressure and energy density from gamma function
#   Author: Alana Gudinas

############################################################################

# Script computes pressure and energy density from gamma coefficients function
# for use in MinimizeErrorPiecewise.

############################################################################

import scipy
import numpy as np
from math import log, exp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.integrate as integrate

import PiecewiseEOS_Read as pp
from NeutronStarEOSlibrary import * 


def polynomial(x,A):
    n = len(A)
    y = A[n-1]*x
    for i in range (n-2,0,-1):
        y = (y + A[i])*x
    return y + A[0]


class CompPresEps(object):
    
    def __init__(self,coefvec):
        self.coefvec = coefvec
    
    def Gamma_function(self,x):
        if x < 0:
            return self.coefvec[0]
        return polynomial(x,self.coefvec)
        
    def PressureFromGammaMR(self,rho):
        if rho < 1.e-20:
            return 1.e-20
        x = log(rho/DensNucCode[0])
        if x < 0:
            return initialPressure*exp(self.coefvec[0]*x)
        Gamma_function_integral = integrate.quadrature(lambda q: self.Gamma_function(q), 0, x, tol = 1e-30, rtol = 1e-3,maxiter=14, vec_func=False, miniter=14)[0]
        return initialPressure*exp(Gamma_function_integral)
    
    def EpsilonFromGammaMR(self,rho):
        if rho < 1.e-20:
            return 1.e-20
        x = log(rho/DensNucCode[0])
        if x < 0:
            return initialEpsilon*exp(x) + initialPressure * ((exp(G0*x)-exp(x))/(G0-1))
        Pressure_function_integral = integrate.quadrature(lambda q: self.PressureFromGammaMR(DensNucCode[0]*exp(q))*exp(-q), 0, x, tol = 1e-30, rtol = 1e-3,maxiter=14, vec_func=False, miniter=14)[0]
        return initialEpsilon*exp(x) + exp(x)*Pressure_function_integral

def main(argv):
    cpe = CompPresEps(testco)

