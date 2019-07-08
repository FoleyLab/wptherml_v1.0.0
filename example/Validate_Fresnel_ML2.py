#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 05:14:09 2019

@author: jay
"""

from wptherml.wpml import multilayer
from matplotlib import pyplot as plt
import numpy as np

''' In the following examples, we will compute the reflectivity, transmissivity, and absorptivity/emissivity
    of simple noble metal films.  Validation data was obtained using S. Byrnes tmm program, 
    which can be obtained from: https://github.com/sbyrnes321/tmm '''

### of the multilayer structure you want to simulate
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        ### temperature of the structure - relevant for all thermal applications
        ### value is stored in attribute self.T
        'Temperature': 300,
        ### gold film with air above and belo
        'Material_List': ['SiO2', 'Au', 'Air'],
        ### gold film is 50 nm thick
        'Thickness_List': [0, 50e-9, 0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [200e-9, 1000e-9, 5000]
        }

gold_film = multilayer(structure)




### use the angular_fresnel method of the multilayer class to compute 
### the reflectivity vs angle for the gold_slab structure at a specified wavelength (in SI units)
gold_film.angular_fresnel(616e-9)

### get the validation data for case 2 - 50 nm Au at 616 nm
gold_film.validation_option = 2
gold_film.get_validation_data()


plt.plot(gold_film.theta_array*180./np.pi, gold_film.r_vs_theta, 'red', label="WPTherml")
plt.plot(gold_film.valid_theta_array, gold_film.valid_ref_vs_theta, 'b--', label="Validation")
plt.xlabel("Angle of Incidence (degrees)")
plt.ylabel("Reflectivity")
plt.legend()
plt.show()

