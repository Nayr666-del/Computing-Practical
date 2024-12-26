# -*- coding: utf-8 -*-
"""
Code to produce Plot 4: Mass-Radius Relation

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

plt.plot(np.array(mass_array)/(1.99*1e30),np.array(r_array)/(696340000),'r-',label = 'relativistic')
plt.plot(np.array(mass_nonrel)/(1.99*1e30),np.array(r_nonrel)/(696340000),'b-',label = 'non-relativistic')
plt.axis([0,3,0,0.07])
y= np.linspace(0,0.07,100)
x = np.ones(100)
plt.plot(mass_array[-1]/(1.99*1e30)*x,y,'g--')
plt.xlabel('Mass/Solar Mass')
plt.ylabel('Radius/Solar Radius')
plt.title('White Dwarf Mass-Radius Relation')
plt.grid()
plt.legend()
print(mass_array[-1]/(1.99*1e30))