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

def Prad_prime(dim, TE_prime_p, TE_prime_s, lam, theta, w):
    dlam = np.abs(lam[0] - lam[1])
    grad = np.zeros(dim)
    ### loop over dof
    for i in range(0,dim):
        x = 0
        ### loop over angle
        for j in range(0,len(w)):
            prad_som = 0.
            ### loop over wavelength
            for k in range(0,len(lam)):
                prad_som = prad_som + (0.5*TE_prime_p[i,j,k] + 0.5*TE_prime_s[i,j,k])*dlam
            x = x + prad_som*(np.sin(theta[j])*w[j])
        grad[i] = 2*np.pi*x
        
    return grad

  
#self.atmospheric_power_val = coolibglib.Patm(self.emissivity_array_p, self.emissivity_array_s, self.T_amb, self.lambda_array, self.t, self.w)

def Patm(EPS_P, EPS_S, T_amb, lam, theta, w):
    
    dlam = np.abs(lam[1] - lam[0])
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

def Patm_prime(dim, eps_prime_p, eps_prime_s, T_amb, lam, theta, w):
    
    dlam = np.abs(lam[1]- lam[0])
    ### get normal atmospheric transmissivity
    atm_T = datalib.ATData(lam)
    ### get BB spectrum associated with ambient temperature
    BBs = datalib.BB(lam, T_amb)
    ### initialize grad
    grad = np.zeros(dim)
    
    for i in range(0,dim):
        x = 0
        for j in range(0,len(w)):
            patm_prime = 0
            angular_mod = 1./np.cos(theta[j])
            for k in range(0,len(lam)):
                patm_prime = patm_prime + 0.5*(eps_prime_p[i,j,k]+eps_prime_s[i,j,k])*BBs[k]*np.cos(theta[j])*(1-atm_T[k]**angular_mod)*dlam
            x = x + patm_prime * np.sin(theta[j]) * w[j]
        grad[i] = 2*np.pi*x
        
    return grad

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

### P sun!
def Psun_prime(grad_list, theta_sun, lam, n, d):
    ### length of arrays 
    ### length of gradient list
    dim = len(grad_list)
    g = np.zeros(dim)
    ### numer of wavelengths
    n_lam = len(lam)
    ### number of layers
    n_layer = len(d)
    
    ### get Am1.5 spectrum
    AM = datalib.AM(lam)
    ### variables to hold emissivity for s- and p-polarized light

    ### allocate array for refractive index at a particular wavelength
    nc = np.zeros(n_layer, dtype=complex)
    
    dl = np.abs(lam[1]-lam[0])
    for i in range(0,n_lam):
        for j in range(0,n_layer):
            nc[j] = n[j][i]
              
                
        k0 = np.pi*2/lam[i]
        ### get p-polarized transfer matrix for this k0, th, pol, nc, and d
        Mp = tmm.tmm_grad(k0, theta_sun, 'p', nc, d, grad_list)
        Ms = tmm.tmm_grad(k0, theta_sun, 's', nc, d, grad_list)
        
        ### p-polarized arrays - normal arrays
        Mp21 = Mp["M21"]
        Mp11 = Mp["M11"]
        ### p-polarized arrays - gradient arrays
        Mp21p = Mp["Mp"][:,1,0]
        Mp11p = Mp["Mp"][:,0,0]
                
        ### s-polarized arrays - normal arrays
        Ms21 = Ms["M21"]
        Ms11 = Ms["M11"]
        ### s-polarized arrays - gradient arrays
        Ms21p = Ms["Mp"][:,1,0]
        Ms11p = Ms["Mp"][:,0,0]
                
        ### p-polarized amplitudes
        rp = Mp21/Mp11
        rp_star = np.conj(rp)
        tp = 1./Mp11
        tp_star = np.conj(tp)
                
        ### s-polarized amplitudes
        rs = Ms21/Ms11
        rs_star = np.conj(rs)
        ts = 1./Ms11
        ts_star = np.conj(ts)
                
        ### get incident/final angle... will be independent of polarization
        ti = Mp["theta_i"]
        tL = Mp["theta_L"]
        fac = nc[len(d)-1]*np.cos(tL)/(nc[0]*np.cos(ti))
        
        ### now we need to loop over the number of elements we are differentiating with respect
        ### to!
        for k in range(0,dim):
            ### p-polarized primed quantities
            rp_prime = (Mp11*Mp21p[k] - Mp21*Mp11p[k])/(Mp11*Mp11)
            tp_prime = -Mp11p[k]/(Mp11*Mp11)
            rp_prime_star = np.conj(rp_prime)
            tp_prime_star = np.conj(tp_prime)
                    
            ### s-polarized primed quantities
            rs_prime = (Ms11*Ms21p[k] - Ms21*Ms11p[k])/(Ms11*Ms11)
            ts_prime = -Ms11p[k]/(Ms11*Ms11)
            rs_prime_star = np.conj(rs_prime)
            ts_prime_star = np.conj(ts_prime)
                    
            ### p-polarized R, T, and epsilon prime
            Rp_prime = rp_prime * rp_star + rp * rp_prime_star
            Tp_prime = (tp_prime * tp_star + tp * tp_prime_star)*fac
            eps_p_prime = 0 - Rp_prime - Tp_prime
            
            ### s-polarized R, T, and epsilon prime
            Rs_prime = rs_prime * rs_star + rs * rs_prime_star
            Ts_prime = (ts_prime * ts_star + ts * ts_prime_star)*fac
            eps_s_prime = 0 - Rs_prime - Ts_prime
            
            g[k] = g[k] + 0.5*(eps_p_prime + eps_s_prime) * AM[i] * dl
            
            
    
    return g


def Pwr_cool(lam):
    Pcool = Prad()-Patm(lam, theta, T)-Psun(lam)
    return Pcool
