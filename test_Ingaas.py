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
pv2 = {
        ### actual materials the structure is made from... note terminal layers are air and
   	### top-side layer (layer upon which light is incident) is SiO2.
        ### Refractive index values are stored in the attribute self.n
        'Material_List': ['Air', 'InGaAs', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 100e-9,  0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [400e-9, 6000e-9, 1000]
        }

pv1 = {
        ### actual materials the structure is made from... note terminal layers are air and
   	### top-side layer (layer upon which light is incident) is SiO2.
        ### Refractive index values are stored in the attribute self.n
        'Material_List': ['Air', 'Si', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 100e-9,  0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [400e-9, 6000e-9, 1000]
        }


em = {
        ### actual materials the structure is made from... note terminal layers are air and
   	### top-side layer (layer upon which light is incident) is SiO2.
        ### Refractive index values are stored in the attribute self.n
        'Material_List': ['Air', 'W', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 900e-9,  0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [400e-9, 6000e-9, 1000],
        'Temperature': 2000
        }


### create the instance called coated_au_film

pv_1 = multilayer(pv1)
pv_2 = multilayer(pv2)
emitter = multilayer(em)
emitter.thermal_emission()

### create a plot of the reflectivity of the coated au film - use red lines
### the wavelengths are stored in SI units so we will multiply by 1e9 to 
### plot them in nanometers
#plt.plot(1e9*coated_au_film.lambda_array, coated_au_film.reflectivity_array, 'red')
#plt.show()

jsc1, jsc2 = stpvlib.jsc_multi(emitter.lambda_array, emitter.thermal_emission_array, pv_1.emissivity_array, pv_2.emissivity_array, pv_1.transmissivity_array )
print(jsc1)
print(jsc2)

