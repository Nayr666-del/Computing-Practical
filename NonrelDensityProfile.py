# -*- coding: utf-8 -*-
"""
Code to produce Plot 1a): Non-relativisitc density profile
"""

import matplotlib.pyplot as plt
import numpy as np
from CO31_utilites import get_density

r = np.linspace(1,5e7,1000)
mass_nonrel, radius_nonrel = get_density(1e12, r, False)

plt.plot(r,mass_nonrel[:,0],'b-',label = 'non-relativistic')
plt.axis([1,0.5e7,0,1e12])
plt.xlabel('Radius/m')
plt.ylabel('Density/(kg/m^3)')
plt.title('Density Profile')
plt.legend()
plt.grid()
