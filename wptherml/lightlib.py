# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 12:32:14 2018

@author: varnerj

This function passed the following tests:
    Predicted luminous efficiency of 14% for a blackbody at 6,600 K in 
    agreement with http://www.ccri.edu/physics/keefe/light.htm
    
    Wikipedia has tabulated values of various standards for Luminous Efficiency and Luminous Efficacy:
    https://en.wikipedia.org/wiki/Luminous_efficacy

"""

###Lightlib
from wptherml.numlib import numlib
from wptherml.datalib import datalib
import numpy as np
from matplotlib import pyplot as plt

    
def luminous_efficiency(lam, thermal_emission_array):
    upper = np.amax(lam)
    ph = datalib.PhLum(lam)
    num = datalib.PhLum(lam)*thermal_emission_array
    numerator = numlib.Integrate(num, lam, 0, upper )
    den = thermal_emission_array
    denominator = numlib.Integrate(den, lam, 0, upper)
    luminous_efficiency_value = (numerator/denominator)
    return luminous_efficiency_value

def luminous_efficacy(lam, thermal_emission_array):
    le = luminous_efficiency_value(lam, thermal_emission_array)
    luminous_efficacy_value = 683*le
    return luminous_efficacy_value

def thermal_emission_ideal_source(lam, T):
    ### compute blackbody spectrum at current T
    rho = datalib.BB(lam, T)
    ### get photopic luminosity function
    ph  = datalib.PhLum(lam)
    ### ideal thermal emission is product of the two
    thermal_emission_ideal_source = ph * rho
    return thermal_emission_ideal_source_array
