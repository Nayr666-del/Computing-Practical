# -*- coding: utf-8 -*-
"""
Code to produce Plot 3b): Radius-Density Relation

Author: Hongyi Xiong jesu4590
"""
import matplotlib.pyplot as plt
import numpy as np
from CO31_utilites import get_density

rho_array = np.logspace(6,16,30)
r_array = []
r_nonrel = []
mass_array = []
mass_nonrel = []
r = np.linspace(1,5e7,10000)
for rho0 in rho_array:
    massdensity, radius= get_density(rho0, r, True)
    massdensity_nonrel, radius_nonrel = get_density(rho0, r, False)
    r_array.append(radius)
    r_nonrel.append(radius_nonrel)
    mass_array.append(massdensity[-1][1])
    mass_nonrel.append(massdensity_nonrel[-1][1])
    

plt.loglog(rho_array,r_nonrel,'bx',label='non-relativistic')
plt.loglog(rho_array,r_array,'r^-',label = 'relativistic')
plt.loglog(rho_array,r_nonrel[0]*(rho_array/rho_array[0])**(-1/6),'g-',label='Power Law')
plt.legend()
plt.title('Radius-Central Density Relation')
plt.xlabel('Central Density/(kg/m^3)')
plt.ylabel('Star Radius/m')
plt.grid()
plt.axis([1e6,1e14,1e5,1e9])