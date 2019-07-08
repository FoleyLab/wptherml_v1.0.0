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

### dictionary that stores basic properties 
### of the multilayer structure you want to simulate
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        ### temperature of the structure - relevant for all thermal applications
        ### value is stored in attribute self.T
        'Temperature': 300,
        ### gold film with air above and belo
        'Material_List': ['Air', 'Au', 'Air'],
        ### gold film is 50 nm thick
        'Thickness_List': [0, 50e-9, 0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [200e-9, 1000e-9, 5000]
        }

gold_film = multilayer(structure)

### get validation data for first validation case - 50 nm Au film
gold_film.validation_option = 1
gold_film.get_validation_data()



plt.plot(1e9*gold_film.lambda_array, gold_film.reflectivity_array, 'red', label='WPTherml')
plt.plot(1e9*gold_film.valid_lambda_array, gold_film.valid_reflectivity_array, 'b--', label='Validation')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectivity")
plt.legend()
plt.show()

### plot wptherml transmissivity with red lines, validation data with blue dashed lines
plt.plot(1e9*gold_film.lambda_array, gold_film.transmissivity_array, 'red', label='WPTherml')
plt.plot(1e9*gold_film.valid_lambda_array, gold_film.valid_transmissivity_array, 'b--', label='Validation')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Transmissivity")
plt.legend()
plt.show()

### plot wptherml emissivity with red lines, validation data with blue dashed lines
plt.plot(1e9*gold_film.lambda_array, gold_film.emissivity_array, 'red', label='WPTherml')
plt.plot(1e9*gold_film.valid_lambda_array, gold_film.valid_emissivity_array, 'b--', label='Validation')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Emissivity")
plt.legend()
plt.show()

