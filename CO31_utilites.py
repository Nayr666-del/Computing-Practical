# -*- coding: utf-8 -*-
"""
A List of utils functions for computing white dwarf structure and producing the plots

Author: Hongyi Xiong jesu4590
"""

import numpy as np

def f_rel(y,r):
# Author: Hongyi Xiong, Date: 18/10/2024
# Computes the relativistic derivatives for the problem
# Input:
# * y: The current state, given by a (1 x M) row vector. In this case it is [density, mass]
# * r: the radius at which the derivative is evaluated
# Output:
# * df: (1 x M) vector containing the derivatives in each element
    mp = 1.672 * 10**(-27) # Proton Mass
    me = 9.11 * 10**(-31) # Electron Mass
    G = 6.673 * 10**(-11) # Gravitation Constant
    h = 6.626 * 10**(-34) # Planck's Constant
    pi = 3.1415926
    c = 299792458
    sinh_theta_F=np.sign(y[0])*h/(2*pi)*(3*pi**2)**(1/3)/(2**(1/3)*me*c*mp**(1/3))*(np.abs(y[0]))**(1/3)
    theta_F = np.arcsinh(sinh_theta_F)
    a = np.cosh(4*theta_F)/(8*np.cosh(theta_F))-np.cosh(2*theta_F)/(2*np.cosh(theta_F))+3/(8*np.cosh(theta_F))
    b = np.sign(y[0])*(9*2**(1/3)*mp**(1/3)*h**2)/(4*me**3*c**4*(3*pi**2)**(1/3))*(np.abs(y[0]))**(2/3)
    df1 = -b/a*G*y[0]*y[1]/r**2
    df2 = 4*pi*r**2*y[0]
    
    return(np.array([df1,df2]))

def f_nonrel(y,r):
# Author: Hongyi Xiong, Date: 18/10/2024
# Computes the derivatives for the problem
# Input:
# * y: The current state, given by a (1 x M) row vector. In this case it is [density, mass]
# * r: the radius at which the derivative is evaluated
# Output:
# * df: (1 x M) vector containing the derivatives in each element
    mp = 1.672 * 10**(-27) # Proton Mass
    me = 9.11 * 10**(-31) # Electron Mass
    G = 6.673 * 10**(-11) # Gravitation Constant
    h = 6.626 * 10**(-34) # Planck Constant
    pi = 3.1415926
    df1 = -48*me*mp**(5/3)*pi**(2/3)/(2**(1/3)*3**(2/3)*h**2)*G*np.abs(y[0])**(1/3)*y[1]/r**2
    df2 = 4*pi*r**2*y[0]
    
    return(np.array([df1,df2]))

def ode_solver_rk(f, y0, r_series):
# Author: Hongyi Xiong, Date: 18/10/2024
# Solve ODE Problem, i.e. dy/dt = f(y,t)
# Solve ODE problem, i.e. dy/dt = f(y,t), using Runge−Kutta algorithm.
# Input:
# ∗ f: a function that receives the current state, y, and the current position/time, t,
# and returns the derivative value of the state, dy/dt.
# ∗ y0: the initial state of the system, given in a row matrix (1 x M).
# ∗ t: vector of position/time steps with length N where the values of y will be returned.
#
# Output:
# ∗ y: (N x M) matrix that contains the values of y at every position/time step
# columns correspond to the position/time and rows to the element of y.
#
# Constraints:
# ∗ Do not use the built−in ODE solvers of Python libraries
    delta_r = r_series[1]-r_series[0]
    y = []
    for r in r_series:
        y.append(list(y0))
        f0 = f(y0,r)
        f1 = f(y0 + delta_r * f0 / 2,r + delta_r / 2)
        f2 = f(y0 + delta_r * f1 / 2,r + delta_r / 2)
        f3 = f(y0 + delta_r * f2,r + delta_r)
        
        y1 = y0 + (f0 + 2*f1 + 2*f2 + f3)*delta_r/6
        
        
        y0 = y1
    
    return y

def get_density(rho0 , r , rel):
# Author: Hongyi Xiong , Date: 18/10/2024
# Obtain the density, rho, as function of the radial distance, r, using the implemented
# ODE solver and the non−relativistic or relativistic equation.
# Input:
# ∗ rho0: the central density at r = 0.
# ∗ r: the grid points of radial distance where the density is calculated in form of
# a vector with N elements.
# ∗ rel: boolean distinguishing relativistic and non−relativistic cases.
 
# Output:
# ∗ rho: an N−element vector that contains the density at the radial grid points given
# in r.
# ∗ mass: the cumulative mass of the white dwarf from r=0 to the radial grid point given
# in r (a vector with N elements).
    mass0 = 4*np.pi/3*rho0
    y0 = np.array([rho0,mass0])
    if rel:
        y_final = np.array(ode_solver_rk(f_rel,y0,r))
    else:
        y_final = np.array(ode_solver_rk(f_nonrel,y0,r))
    x = 0
    for i in range(len(r)):
        if y_final[i][0]<=0:
            x = i
            break
    for i in range(x,len(r)):
        y_final[i][0] = 0.0
        y_final[i][1] = y_final[x-1][1]
    radius = r[x-1]    
    return y_final, radius
