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




def E_atm(theta, lam):
    ### emissivity starts with 1 - T of atmosphere 
    T = datalib.ATData(lam)
    ### angular part is the emissivity raised to 1/cos(theta):
    b = 1./np.cos(theta)
    Eatm = (1-T)**b
    return Eatm




def Prad(TEP, TES, lam, theta, w):
    dlam = np.abs(lam[0] - lam[1])
    x = 0
    for i in range(0,len(w)):
        prad_som = 0
        for j in range(0,len(lam)):
            prad_som = prad_som + (0.5*TEP[i][j] + 0.5*TES[i][j])*dlam
        x = x +prad_som*(np.sin(theta[i])*w[i])
    return 2*np.pi*x

  
#self.atmospheric_power_val = coolibglib.Patm(self.emissivity_array_p, self.emissivity_array_s, self.T_amb, self.lambda_array, self.t, self.w)

def Patm(EPS_P, EPS_S, T_amb, lam, theta, w):
    
    dlam = np.abs(lam[0] - lam[1])
    ### Get normal atmospheric transmissivity
    atm_T = datalib.ATData(lam)
    ### Get BB spectrum associated with ambient temperature
    BBs = datalib.BB(lam, T_amb)
    
    x = 0
    for i in range(0,len(w)):
        patm_som = 0
        angular_mod = 1./np.cos(theta[i])
        for j in range(0,len(lam)):
            patm_som = patm_som + (0.5*EPS_P[i][j] + 0.5*EPS_S[i][j]) * BBs[j] * np.cos(theta[i]) * (1 - atm_T[j]**angular_mod) * dlam
        x = x + patm_som * np.sin(theta[i]) * w[i]
    return 2*np.pi*x

### P sun!
def Psun(theta_sun, lam, n, d):
    ### length of arrays 
    n_lam = len(lam)
    n_layer = len(d)
    
    ### get Am1.5 spectrum
    AM = datalib.AM(lam)
    ### variables to hold emissivity for s- and p-polarized light
    emissivity_s = 0.
    emissivity_p = 0.
    ### allocate array for refractive index at a particular wavelength
    nc = np.zeros(n_layer, dtype=complex)
    
    ### compute emissivity and integrate all at once
    P_sun_sum = 0.
    dl = np.abs(lam[1]-lam[0])
    for i in range(0,n_lam):
        for j in range(0,n_layer):
            nc[j] = n[j][i]
                
        k0 = np.pi*2/lam[i]
        ### get p-polarized transfer matrix for this k0, th, pol, nc, and d
        Mp = tmm.tmm(k0, theta_sun, 'p', nc, d)
        ### get s-polarized transfer matrix for this k0, th, pol, nc, and d
        Ms = tmm.tmm(k0, theta_sun, 's', nc, d)
        ### get t amplitudes
        tp = 1./Mp["M11"]
        ts = 1./Ms["M11"]
            
        ### get incident/final angles
        tpi = Mp["theta_i"]
        tpL = Mp["theta_L"]
        tsi = Ms["theta_i"]
        tsL = Ms["theta_L"]
        ### get geometric factor associated with transmission
        facp = nc[n_layer-1]*np.cos(tpL)/(nc[0]*np.cos(tpi))
        facs = nc[n_layer-1]*np.cos(tsL)/(nc[0]*np.cos(tsi))
        ### get reflection amplitude
        rp = Mp["M21"]/Mp["M11"]
        rs = Ms["M21"]/Ms["M11"]
        ### get Reflectivity
        Rp = np.real(rp * np.conj(rp))
        Rs = np.real(rs * np.conj(rs))
        ### get Transmissivity
        Tp = np.real(tp * np.conj(tp) * facp)
        Ts = np.real(ts * np.conj(ts) * facs)
        emissivity_p = 1 - Rp - Tp
        emissivity_s = 1 - Rs - Ts
        
        P_sun_sum = P_sun_sum + 0.5*(emissivity_p + emissivity_s) * AM[i] * dl
    
    return P_sun_sum



def Pwr_cool(lam):
    Pcool = Prad()-Patm(lam, theta, T)-Psun(lam)
    return Pcool
