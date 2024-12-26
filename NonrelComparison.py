# -*- coding: utf-8 -*-
"""
Code to create Plot 1b): Comparison of different central densities

Author: Hongyi Xiong jesu4590
"""
import matplotlib.pyplot as plt
import numpy as np
from CO31_utilites import get_density

r = np.linspace(1,5e7,1000)
mass_nonrel1, radius_nonrel = get_density(1e6, r, False)
mass_nonrel2, radius_nonrel = get_density(1e8, r, False)
mass_nonrel3, radius_nonrel = get_density(1e10, r, False)
mass_nonrel4, radius_nonrel = get_density(1e12, r, False)
mass_nonrel5, radius_nonrel = get_density(1e14, r, False)

plt.plot(r,mass_nonrel1[:,0],'r-',label = '1e6')
plt.plot(r,mass_nonrel2[:,0],'b-',label = '1e8')
plt.plot(r,mass_nonrel3[:,0],'g-',label = '1e10')
plt.plot(r,mass_nonrel4[:,0],'c-',label = '1e12')
plt.plot(r,mass_nonrel5[:,0],'m-',label = '1e14')
plt.axis([1,5e7,1,1e15])
plt.yscale('log')
plt.xlabel('Radius/m')
plt.ylabel('Density/(kg/m^3)')
plt.title('Density Profiles for differing Central Density')
plt.legend(loc=1)
plt.grid()




