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

def Lum_efficiency(lam, TE):
    upper = np.amax(lam)
    ph = datalib.PhLum(lam)
    num = ph*TE
    numerator = numlib.Integrate(num, lam, 0, upper )
    den = TE
    denominator = numlib.Integrate(den, lam, 0, upper)
    return (numerator/denominator)

def normalized_power(lam, TE, BB):
    upper = np.amax(lam)
    ph = datalib.PhLum(lam)
    num = ph*TE
    den = ph*BB
    numerator = numlib.Integrate(num, lam, 0, upper)
    denominator = numlib.Integrate(den, lam, 0, upper)
    return numerator/denominator

def lum_efficiency_filter(lam, BBs, emissivity, transmissivity):
    upper = np.amax(lam)
    ph = datalib.PhLum(lam)
    ### total observed thermal emission is BB spectrum * emissivity of emitter * transmissivity of filter
    TE = BBs * emissivity * transmissivity
    num = ph * TE
    numerator = numlib.Integrate(num, lam, 0, upper)
    denominator = numlib.Integrate(TE, lam, 0, upper)
    return (numerator/denominator)

def lum_efficiency_filter_prime(dim, lam, BBs, emissivity, transmissivity, transmissivity_prime):
        
    ### allocate gradient vector
    grad = np.zeros(dim)
    ### get data that will not change for each element of the gradient first!
    upper = np.amax(lam)
    ph = datalib.PhLum(lam)
    TE = BBs * emissivity * transmissivity
    TE_prime = np.zeros(len(lam))
    num = ph*TE
    numerator = numlib.Integrate(num, lam, 0, upper )
    denominator = numlib.Integrate(TE, lam, 0, upper)
    
    ### now loop through elements of gradient_list and fill in elements of grad
    for i in range(0,dim):
        TE_prime = BBs * emissivity * transmissivity_prime[i,:]
        num_prime = ph * TE_prime
        numerator_prime = numlib.Integrate(num_prime, lam, 0, upper)
        denominator_prime = numlib.Integrate(TE_prime, lam, 0, upper)
        grad[i] = (denominator*numerator_prime - numerator*denominator_prime)/(denominator**2)
    
    return grad

def lum_efficiency_prime(dim, lam, emissivity, emissivity_prime, BBs):
    
    ### allocate gradient vector
    grad = np.zeros(dim)
    ### get data that will not change for each element of the gradient first!
    upper = np.amax(lam)
    ph = datalib.PhLum(lam)
    TE = emissivity*BBs
    TE_prime = np.zeros(len(lam))
    num = ph*TE
    numerator = numlib.Integrate(num, lam, 0, upper )
    denominator = numlib.Integrate(TE, lam, 0, upper)
    
    ### now loop through elements of gradient_list and fill in elements of grad
    for i in range(0,dim):
        TE_prime = BBs*emissivity_prime[i,:]
        num_prime = ph*TE_prime
        numerator_prime = numlib.Integrate(num_prime, lam, 0, upper)
        denominator_prime = numlib.Integrate(TE_prime, lam, 0, upper)
        grad[i] = (denominator*numerator_prime - numerator*denominator_prime)/(denominator**2)
    
    return grad


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


