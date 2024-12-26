# -*- coding: utf-8 -*-
"""
Code to Produce Plot 3a): Mass-density relations

@author: Ryan Xiong
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

plt.loglog(rho_array,mass_array,'r^-',label = 'relativistic')
plt.loglog(rho_array,mass_nonrel,'bx',label = 'non-relativistic')
plt.loglog(rho_array,2.765*1e30*np.ones(len(rho_array)),'c--',label = 'Mass Limit')
plt.loglog(rho_array,mass_nonrel[0]*(rho_array/rho_array[0])**0.5,'g-',label='Power Law')
plt.axis([1e6,1e14,1e28,1e33])
plt.xlabel('Central Density/(kg/m^3)')
plt.ylabel('White Dwarf Mass/kg')
plt.title('Mass-Central Density Relation')
plt.legend()
plt.grid()