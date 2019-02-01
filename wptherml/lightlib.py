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

    
def Lum_efficiency(lam, TE):
    upper = np.amax(lam)
    ph = datalib.PhLum(lam)
    num = datalib.PhLum(lam)*TE
    numerator = numlib.Integrate(num, lam, 0, upper )
    den = TE
    denominator = numlib.Integrate(den, lam, 0, upper)
    return (numerator/denominator)

def Lum_efficacy(lam, TE):
    le = Lum_efficiency(lam, TE)
    efficacy = 683*le
    return efficacy

def IdealSource(lam, T):
    ### compute blackbody spectrum at current T
    rho = datalib.BB(lam, T)
    ### get photopic luminosity function
    ph  = datalib.PhLum(lam)
    ### ideal thermal emission is product of the two
    TE_ideal = ph * rho
    return TE_ideal