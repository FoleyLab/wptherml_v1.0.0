# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 13:36:55 2019

@author: varnerj
"""

from wptherml.wpml import multilayer
from matplotlib import pyplot as plt

### dictionary that stores basic properties 
### of the multilayer structure you want to simulate
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        ### temperature of the structure - relevant for all thermal applications
        ### value is stored in attribute self.T
        'Temperature': 1700,
        ### actual materials the structure is made from
        ### values are stored in the attribute self.n
        'Material_List': ['Air', 'W', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 400e-9, 0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [300e-9, 6000e-9, 1000],
         ### The folloing entry will tell the computer to use the lightbulb functions to
         ### compute properties needed to characterize an incandescent source
         'STPV_EMIT': 1
     
        }

w_ = multilayer(structure)

structure["Material_List"] = ['Air','SiO2','TiO2','SiO2','Al2O3','W','Air'] 
structure["Thickness_List"] = [0,255e-9, 150e-9,255e-9,10e-9,900e-9,0]
w_bragg = multilayer(structure)

structure["Material_List"] = ['Air','W_Al2O3_Alloy','SiO2','TiO2','SiO2','Al2O3','W','Air']
structure["Thickness_List"] = [0,20e-9,255e-9, 150e-9,255e-9,10e-9,900e-9,0]
w_bragg_alloy = multilayer(structure)

### Plot Thermal Emission of bare W slab with a red line
plt.plot(w_.lambda_array*1e9, w_.thermal_emission_array, 'red')
### Plot Thermal Emission of coated W with a blue line
plt.plot(w_bragg.lambda_array*1e9, w_bragg.thermal_emission_array, 'blue')
###
plt.plot(w_bragg_alloy.lambda_array*1e9, w_bragg_alloy.thermal_emission_array, 'purple')
### Plot blackbody spectrum with a black line
plt.plot(w_.lambda_array*1e9, w_.BBs, 'black')
plt.show()





