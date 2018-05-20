############################################################################

#   Equation of state computations from Lindblom paper
#   Author: Alana Gudinas

############################################################################

# Defines functions to compute equation of state as per Lindblom paper.

############################################################################

import scipy
import numpy as np
from math import log
from math import exp
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from numpy import linspace

from NeutronStarEOSlibrary import *


def gamma(p):
    y = log(p/p0)
    return exp(g0 + g1*y + g2*(y**2)+ g3*(y**3))

def gamma1(p):
    return 1/gamma(p)

def mu(p):
    mid = integrate.quad(lambda q: (1/q)*gamma1(q), p0, p)[0]
    return exp(-mid)

def combEpsGam(p):
    return mu(p)/gamma(p)

def epsilon(p):
    eps1 = integrate.quad(lambda q: combEpsGam(q), p0, p)[0]
    # epsilon(p) = (eps0/mu(p)) + (1/mu(p))*eps1 
    return (eps0/mu(p)) + (1/mu(p))*eps1 

def main():
    #result = epsilon(1)
    #print(result)
    print(gamma(p0))

x = linspace(1,100,100)
for j in range(0,100):
    x[j] = exp(10.*j/100.)
x = x*p0
epsilon_lind = np.zeros(100)

for j in range(0,100):
    epsilon_lind[j] = (epsilon(x[j]))
    epsilon_lind[j] = log(epsilon_lind[j])
    x[j] = log(x[j])
    
plt.plot(x,epsilon_lind)
plt.show()


