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


'''
def E_atm(T,theta):
    ### emissivity starts with 1 - T of atmosphere 
    a = 1-T
    ### angular part is the emissivity raised to 1/cos(theta)
    b = 1/np.cos(theta)
    Eatm = a**b
    
    return Eatm
'''

'finish Prad function to remove error'

def Prad(TEP, TES, lam, theta, w):
    dlam = np.abs(lam[0] - lam[1])
    x = 0
    for i in range(0,len(w)):
        prad_som = 0
        for j in range(0,len(lam)):
            prad_som = prad_som + ( 0.5*TEP[i][j] + 0.5*TES[i][j])*dlam
        x = x +prad_som*(np.sin(theta[i])*w[i])
    return x

'''   
def Patm(lam, theta, T):
    dlam = np.abs(lam[0] - lam[1])
    x = 0
    for i in range(0,len(w)):
        patm_som = 0
        for j in range(0,len(lam)):
            patm_som = patm_som + ()
 'finish Psun funtion'   
def Psun(A, lam, TES, TEP,theta):
    upper = np.amax(lam)
    AM = datalib.AM(lam)
    E = tmm.abs()
    Psun = numlib.Integrate(AM*E,lam,100e-9, upper)
    return Psun




def Pwr_cool(lam):
    Pcool = Prad()-Patm(lam, theta, T)-Psun(lam)
    return Pcool
   ''' 
