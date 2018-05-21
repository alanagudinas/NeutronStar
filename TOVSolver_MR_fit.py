############################################################################

#   TOV solver: adapted for mass/radius computations
#   Author: Francois Foucart

############################################################################


from numpy import sqrt,pi
from scipy.optimize import brentq
import numpy as np
from math import log, exp
import matplotlib.pyplot as plt

from PiecewiseEOS_Read import *

from PresEpsCompMR import CompPresEps


testrho = 0.002
rhoMax = 0.003


rhostart = 0.0003


def runcomp(coefvec):
    comp = CompPresEps(coefvec)
    return comp

###### Calls equation of state function ####

def EpsilonFromDensity(rho, cpe):
    return cpe.EpsilonFromGammaMR(rho)
    
def PressureFromDensity(rho, cpe):
    return cpe.PressureFromGammaMR(rho)

def ZeroFuncForPressure(rho,P,cpe):
    if(rho<0):
        return 0.
    return PressureFromDensity(rho, cpe)-P


# print EpsilonFromDensity(testrho, comp)

####### Solve for density given pressure ####
def DensityFromPressure(P,cpe):
    if(P<0):
        return 0.
    root = brentq(ZeroFuncForPressure,1.e-16,rhoMax,args=(P,cpe),xtol=1.e-15,rtol=1.e-15)
    return root

###### Computes derivatives of ODE - should not need any change #####
# Assumes y[0]=r_circ, y[1]=m, y[2]=phi, y[3]=P, y[4]=rest mass
    
def TOVderivs(t,y, cpe):
    # Give physical names to variables
    r_p=t
    rbar_p=y[0]
    m_p=y[1]
    # Find mass and energy density
    P_p = y[3]
    rho_p = DensityFromPressure(P_p,cpe)
    eps_p = EpsilonFromDensity(rho_p,cpe)#rho_p*(1.+EpsilonFromDensity(rho_p))
    #Computation of derivatives
    A=1.
    if(r_p>1.e-16):
        A = (1.-2.*m_p/r_p)
    sqA=sqrt(A)
    dydx = 1.*y
    dydx[0] = 1.
    if(r_p>1.e-16):
        dydx[0] = 1./(sqA*(r_p/rbar_p))
    dydx[1]=4.*pi*(r_p**2)*eps_p
    dydx[2]=0.
    dydx[3]=0.
    dydx[4]=0.
    if(m_p>1.e-16):
        dydx[2]=(m_p+4.*pi*P_p*r_p**3)/((r_p**2)*A)
        dydx[3]=-(eps_p+P_p)*dydx[2]
        dydx[4]=4*pi*(r_p**2)*rho_p/sqA
    return dydx

from scipy.integrate import ode

def IntegrateTOVSolution(rho_c,dR,cpe):
    ### Initialize ODE solver ####
    ode_int = ode(TOVderivs)
    #### Initial conditions: central density and pressure ###
    P_c = PressureFromDensity(rho_c,cpe)
    y0 = [0,0.,0.,P_c,0.]
    ode_int.set_initial_value(y0,0).set_f_params(cpe)
    #### Integrate equations ####
    step = 0
    while ode_int.y[3]>0.:
        ode_int.integrate(ode_int.t+dR)
        #if(step%100 == 0):
        # 3print(ode_int.y)
        step = step+1
    res = [ode_int.y[0],ode_int.y[1],ode_int.y[2],ode_int.y[3],ode_int.y[4],ode_int.t]
    return res

import numpy as np
from numpy import zeros as zeros

def testcalc(cpe):
    rho_c = .002
    dR = 1.e-3
    sol = zeros(5)
    sol = IntegrateTOVSolution(rho_c,dR,cpe)
    #### Get output ####
    # Outer radius
    r_f = sol[5]
    # Total gravitational mass
    m_f = sol[1]
    # Outer isotropic radius (different coordinates...)
    r_iso = 0.5*(r_f-m_f+r_f*sqrt(1.-2.*m_f/r_f))
    # Pressure on surface
    P_f = sol[3]
    # Baryon mass (> grav. mass...)
    RestMass = sol[4]
    ##### Output #####
    print("ADM Mass = %g" % m_f)
    print("Radius (isotropic) = %g" % r_iso)
    print("Radius (circular) = %g" % r_f)
    print("Pressure = %g" % P_f)
    print("RestMass = %g" % RestMass)


# originally started with rho0 = .00023330409
def MRfitvec(cpe):
    rho_c = np.arange(0.0003,3.0e-3,0.0002) ## start with smaller vector to test
    n = len(rho_c)
    radius_t = zeros(n)
    mass_t = zeros(n)
    dR = 1.e-3
    for i in range (0,n):
        centralrho = rho_c[i]
        #print("Central density = %g" % rho_c[i])
        sol = zeros(5)
        sol = IntegrateTOVSolution(centralrho,dR,cpe)
        # Outer radius
        r_f = sol[5]
        # Total gravitational mass
        m_f = sol[1]
        # Outer isotropic radius (different coordinates...)
        r_iso = 0.5*(r_f-m_f+r_f*sqrt(1.-2.*m_f/r_f))
        # Pressure on surface
        P_f = sol[3]
        # Baryon mass (> grav. mass...)
        RestMass = sol[4]
        radius_t[i] = sol[5]
        #print radius_t[i]
        mass_t[i] = sol[1]
        #print mass_t[i]
    return mass_t,radius_t

def main(coefvec):
    global comp
    comp = runcomp(coefvec)
    mass,radius = MRfitvec(comp)
    plt.plot(radius,mass)
    plt.show()
    return mass, radius



