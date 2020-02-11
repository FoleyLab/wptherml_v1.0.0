# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 13:57:44 2020

@author: foleyj10
"""

from wptherml.wpml import multilayer
from matplotlib import pyplot as plt
from wptherml import stpvlib
### dictionary that stores basic properties 
### of the multilayer structure you want to simulate
structure = {
        ### actual materials the structure is made from... note terminal layers are air and
   	### top-side layer (layer upon which light is incident) is SiO2.
        ### Refractive index values are stored in the attribute self.n
        'Material_List': ['Air', 'SiO2', 'TiO2', 'InGaAs', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 100e-9, 50e-9, 20e-9,  0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [400e-9, 6000e-9, 1000]
        }

### create the instance called coated_au_film
coated_au_film = multilayer(structure)

### create a plot of the reflectivity of the coated au film - use red lines
### the wavelengths are stored in SI units so we will multiply by 1e9 to 
### plot them in nanometers
#plt.plot(1e9*coated_au_film.lambda_array, coated_au_film.reflectivity_array, 'red')
#plt.show()

stpvlib.jsc_multi(coated_au_film.lambda_array)

