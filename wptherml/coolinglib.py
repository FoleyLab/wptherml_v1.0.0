# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 13:14:46 2018

@author: varnerj
"""

###cooling lib
from wptherml import tmm
from wptherml.datalib import datalib
from wptherml import stpvlib
import numpy as np
import matplotlib.pyplot as plt
from wptherml.numlib import numlib

''' 
   
    these lines are outside of all functions
    and probably should be deleted

Ep= tmm.abs(k0, theta0, P, nA, tA)
Es=tmm.abs(k0, theta0, s, nA, tA)

'''

def Eatm(T,theta):
    a = 1-T
    b = 1/np.cos(theta)
    Eatm = a**b
    return Eatm

def Prad():
    
def Patm(lam, theta, T):
    upper = np.amax(lam)
    BB = stpvlib.BB(lam, T)
    A = TMM.abs()
    TE = A*BB
    trig = np.sin(theta)*np.cos(theta)
    ET = (Ep+Es)/2
    Patm=np.pi*2*numlib*Integrate(trig,theta,0, np.pi/2)*Integrate(TE,lam,100e-9,upper)*ET
    return Patm
    
def Psun(lam):
    upper = np.amax(lam)
    AM = datalib.AM(lam)
    E = TMM.abs()
    numlib.Integrate(AM*E,lam,100e-9, upper)




def Pwr_cool(lam):
    Pcool = Prad()-Patm(lam, theta, T)-Psun(lam)
    return Pcool
    
