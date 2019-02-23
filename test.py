#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:35:19 2019

@author: jay
"""

from wptherml.wpml import multilayer
from matplotlib import pyplot as plt
from wptherml.datalib import datalib

### dictionary that stores basic properties 
### of the multilayer structure you want to simulate
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        ### temperature of the structure - relevant for all thermal applications
        ### value is stored in attribute self.T
        'Temperature': 300,
        ### actual materials the structure is made from
        ### values are stored in the attribute self.n
        'Material_List': ['Air', 'SiO2', 'HfO2', 'SiO2', 'HfO2', 'SiO2', 'HfO2', 'SiO2', 'Ag', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 230e-9, 485e-9, 688e-9, 13e-9, 73e-9, 34e-9, 54e-9, 200e-9, 0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [300e-9, 20000e-9, 5000],
        
        'EXPLICIT_ANGLE': 1,
         ### The folloing entry will tell the computer to use the lightbulb functions to
         ### compute properties needed to characterize an incandescent source
         'LIGHTBULB': 1
     
        }

### create the instance called glass_slab
### Create figure that resembles Figure 2 a in Fan et al passive cooling paper:  
w_slab = multilayer(structure)
AM = datalib.AM(w_slab.lambda_array)

plt.plot(w_slab.lambda_array*1e6, w_slab.emissivity_array, 'blue')
plt.plot(w_slab.lambda_array*1e6, AM/(1.4*1e9), 'red')
plt.xlim(0.3,2.5)
plt.show()

### Create figure that resembles Figure 2b in Fan et al passive cooling paper: 
#T_atm = datalib.ATData(w_slab.lambda_array)
#plt.plot(w_slab.lambda_array*1e6, w_slab.emissivity_array, 'red')
#plt.plot(w_slab.lambda_array*1e6, T_atm, 'blue')
#plt.xlim(2.5,20)
#plt.show()



w_slab.cooling_power()

print(w_slab.radiative_power_val)
print(w_slab.solar_power_val)
print(w_slab.atmospheric_power_val)
print(w_slab.cooling_power_val)
