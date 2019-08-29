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


### This figure of merit takes the overlap of the coupled J-Agg / ML structure
### with a 15 nm J-Agg layer and normalizes it to the overlap of the 15 nm J-Agg absorbance with itself
### hence: this is a dimensionless quantity or "enhancement" factor.  If the
### absorbance of the coupled structure is identical to the J-Agg absorbance, the fom will be 1 (no enhancement).
### If the composite absorbance is diminished relative to the lone J-Agg, the fom will be < 1,
### and if it is enhanced relative to lone J-Agg, it will be > 1
def Jagg_enhancement(emissivity, lam):
    upper = np.amax(lam)
    jg = datalib.JAgg_Abs(lam)
    num = jg*emissivity
    numerator = numlib.Integrate(num, lam, 0, upper )
    den = jg*jg
    denominator = numlib.Integrate(den, lam,0, upper)
    return (numerator/denominator)


### Only the numerator has terms that depend on thickness (lone j-agg absorbance is constant with thickness of ML layers)
### so we only need emissivity_prime, there are complicated terms from the quotient of two quantities that depend on thickness
def Jagg_enhancement_prime(dim, lam, emissivity, emissivity_prime):
    grad = np.zeros(dim)
    upper = np.amax(lam)
    #upper = np.amax(lam)
    ### we want the absorption spectrum of the lone J-Agg layer, not just the refractive index
    jg = datalib.JAgg_Abs(lam)
    denominator = jg*jg
    ### denominator intergral only needs to be computed once
    den = numlib.Integrate(denominator, lam, 0, upper)
    for i in range(0,dim):
        numerator = emissivity_prime[i,:] * jg
        num_prime = numlib.Integrate(numerator, lam, 0, upper)
        grad[i] = num_prime / den

    return grad


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


