#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 21:07:22 2018

@author: jay
"""
from wptherml import tmm
from wptherml import colorlib
from wptherml import coolinglib
from wptherml import stpvlib
from wptherml import lightlib
from wptherml.numlib import numlib
from wptherml.datalib import datalib
from matplotlib import pyplot as plt
from scipy import integrate
import numpy as np
#from pkg_resources import resource_string
#import os
class multilayer:
    
    ### initializer
    #def __init__(self, mode, inputfile):
    def __init__(self, args):
        ### set up some default attributes
        #contents = resource_string('wptherml','datalib/AM_15.txt')
        #print(contents)
        #print(os.path.realpath(__file__))
        self.result = 1
        self.mode = args.get('mode')
        #print(" Mode is ",self.mode)
        
        ### There might not always be an input file..
        ### need to make sure we think about to handle this
        self.inputfile = args.get('file')
        
        ### Set default values for all attributes
        ###that are required for actually running
        ### different calculations... will check for
        ### user input on these later
        self.pol = 'p'
        ### default incident angle
        self.theta = 0
        ### default solar angle for coolinglib calculations
        self.theta_sun = 30 * np.pi/180
        ### T_ml is the temperature of the multilayer being modeledd
        self.T_ml = 300
        ### T_amb is the ambient temperature
        self.T_amb = 300
        ### default bandgap wavelength is 2254e-9 m, good value for InGaAsSb
        self.lbg = 2254e-9
        ### default PV is InGaAsSb
        self.PV = "InGaAsSb"
        ### default Temperature of PV cell is 25 deg C 
        ### or 298 K
        self.T_cell = 298
        ### default solar concentration for STPV absorber
        self.solarconc = 600
        ### default is to not use explicit angle dependence of emissivity, etc
        self.explicit_angle = 0
        ### by default, degree of G-L polynomial will not be needed since explicit angle dependence
        ### is not the default but 
        ### we will set it at 7 anyway for now
        self.deg = 7
        
        ### relates to SPP and PA resonances
        self.SPP_Resonance = 0+0j
        self.PA_Resonance = 0+0j
        ### attributes that relate to which quantities should
        ### be computed... the default will be 
        ### only to compute the Fresnel reflection, transmission,
        ### and absorption/emissivity for a structure
        ### all further options must be specified by user
        ### typically stpv emitters will be designed
        ### using different criteria than absorbers
        ### so the two applications will be treated independently
        self.stpv_emitter_calc = 0
        self.stpv_absorber_calc = 0
        ### note there could be a third section for an integrated absorber/emitter
        ### but some thought is required 
        self.cooling_calc = 0
        self.lightbulb_calc = 0
        self.spp_calc = 0
        self.color_calc = 0
        self.fresnel_calc = 1
        self.explicit_angle = 0
        self.resonance = 0
        self.reflective_rgb = np.zeros(3)
        self.thermal_rgb = np.zeros(3)
        self.color_name = 'None'
        
        ### current version only inline_structure method supported
        ### more modes of operation will come in later versions!
        self.inline_structure(args)
            
        ### Now that structure is defined and we have the lambda array, 
        ### allocate other arrays!
        ### Always need normal arrays
        self.reflectivity_array = np.zeros(len(self.lambda_array))
        self.transmissivity_array = np.zeros(len(self.lambda_array))
        self.emissivity_array = np.zeros(len(self.lambda_array))
        self.thermal_emission_array = np.zeros(len(self.lambda_array))
        
        ### if we run validation tests, the following will be used!
        self.valid_lambda_array = []
        self.valid_reflectivity_array = []
        self.valid_transmissivity_array = []
        self.valid_emissivity_array = []
        self.vald_theta_array = []
        self.valid_ref_vs_theta = []
        self.valid_trans_vs_theta = []
        self.valid_emiss_vs_theta = []
        self.validation_option = 1
        
        ### derivative quantities
        self.reflectivity_prime_array = np.zeros(len(self.lambda_array))
        self.emissivity_prime_array = np.zeros(len(self.lambda_array))
        
        
        ### In some cases the user may wish to compute
        ### R, T, or eps vs angle at a specific wavelength
        ### we will allocate three arrays for these cases
        ### with a resolution of 0.5 degrees... i.e. they are small!
        self.r_vs_theta = np.zeros(180)
        self.t_vs_theta = np.zeros(180)
        self.eps_vs_theta = np.zeros(180)
        self.theta_array = np.linspace(0,89.5*np.pi/180, 180)
        ### if we run validation tests, the following will be used!

        ### If users selects explicit_angle option, we 
        ### need arrays for R, T, and eps as a function of angle
        ### and polarization, as well
        if (self.explicit_angle):
            ### range is 0 to thetaC
            a = 0
            b = np.pi/2.
            self.x, self.w = np.polynomial.legendre.leggauss(self.deg)
            self.t = 0.5*(self.x + 1)*(b - a) + a
            self.w = self.w * 0.5 * (b-a)
            
            self.reflectivity_array_p = np.zeros((self.deg,len(self.lambda_array)))
            self.reflectivity_array_s = np.zeros((self.deg,len(self.lambda_array)))
            self.transmissivity_array_p = np.zeros((self.deg,len(self.lambda_array)))
            self.transmissivity_array_s = np.zeros((self.deg,len(self.lambda_array)))
            self.emissivity_array_p = np.zeros((self.deg,len(self.lambda_array)))
            self.emissivity_array_s = np.zeros((self.deg,len(self.lambda_array)))
            self.thermal_emission_array_p = np.zeros((self.deg,len(self.lambda_array)))
            self.thermal_emission_array_s = np.zeros((self.deg,len(self.lambda_array)))


            
        ### Get all far-field Fresnel-related quantities:
        ### Reflectivity, Transmissivity, and Absorptivity/Emissivity spectra
        ### Always call the normal version
        self.fresnel()
        
        ### If user selected explict_angle option, call the EA methods as well
        if (self.explicit_angle):
            ### get the Reflectivity, Transmissivity, and Absorptivity/Emissivity
            ### at the angles from the Gauss-Legendre grid
            self.fresnel_ea()
            ### Get the thermal emission at the angles from the
            ### Gauss-Legendre grid
            self.thermal_emission_ea()
        


        ### stpv_calc cooling_calc lightbulb_calc color_calc all
        ### require BB spectrum / thermal emission spectrum
        if (self.stpv_emitter_calc or self.stpv_absorber_calc or self.cooling_calc or self.lightbulb_calc or self.color_calc):
            ### The ThermalEmission() method automatically calculates the BB spectrum
            ### and stores it to self.BBs
            self.thermal_emission()

        #self.ThermalColor()
        ### now that default quantitites have been calculated, start
        ### looking at optional quantities
        ### want to compute stpv quantities with explicit angle dependence?
        if (self.stpv_emitter_calc and self.explicit_angle):
            
            self.stpv_se_ea()
            self.stpv_pd_ea()
            self.stpv_etatpv_ea()
        ### want no explicit angle dependence?
        elif (self.stpv_emitter_calc):
            
            self.stpv_se()
            self.stpv_pd()
            self.stpv_etatpv()
            
        if (self.stpv_absorber_calc):
        
            if (self.explicit_angle):
                self.stpv_etaabs_ea()
            else:
                self.stpv_etaabs()
            
        if (self.color_calc):

            self.ambient_color()
            self.thermal_color()
            
        if (self.lightbulb_calc):
            
            ### Luminous efficiency and efficacy calcs here
            self.luminous_efficiency()
            #plt.plot(self.lambda_array*1e9, self.thermal_emission_array, 'red', label = 'Lighbulb Emission')
            #plt.plot(self.lambda_array*1e9, self.BBs, 'black', label = 'Blackbody spectrum')
            #plt.xlabel('Wavelength (nm)')
            #plt.ylabel('Spectral Irradiance (W / m^2 / nm / sr)')
            #plt.show()
            
            ### need to validate method for computing luminous efficacy
            #self.luminous_efficacy()
        
        if (self.cooling_calc):
            
            ### will compute the following quantites:
            ### self.radiative_power_val -> at current T_ml, power emitted / unit area
            ### self.solar_power_val -> power absorbed from sun / unit area
            ### self.atmospherical_power_val -> power absorbed from atm at T_amb / unit area
            ### self.cooling_power -> net power flux / unit area from balance of all the above
            self.cooling_power()
            
            AM = datalib.AM(self.lambda_array)
            T_atm = datalib.ATData(self.lambda_array)
            print(self.cooling_power_val," W/m^2 (Total Cooling Power)")
            print(self.radiative_power_val," W/m^2 ((Cooling) Power radiated by structure at ",self.T_ml, "K)")
            print(self.solar_power_val," W/m^2 ((Warming) Power absorbed from sun)")
            print(self.atmospheric_power_val," W/m^2 ((Warming) Power absorbed from atmospheric radiation at ",self.T_amb, "K)")
            plt.xlim(0.3, 2.5)
            plt.plot(self.lambda_array*1e6, self.emissivity_array, 'blue', label = 'Emissivity')
            plt.plot(self.lambda_array*1e6, AM/(1.4*1e9), 'red', label = 'Solar Spectrum')
            plt.xlabel('Wavelength (microns)')
            plt.ylabel('Arb. Units')
            plt.show()
            plt.xlim(2.5, 20)
            plt.plot(self.lambda_array*1e6, T_atm, 'cyan', label = 'Atmospheric Transmissivity')
            plt.plot(self.lambda_array*1e6, self.emissivity_array, 'red', label = 'Emissivity')
            plt.xlabel('Wavelength (microns)')
            plt.ylabel('Arb. Units')
            plt.show()


    def get_validation_data(self):
        Valid_Dict = datalib.read_validation_data(self.validation_option)
        if (self.validation_option==1 or self.validation_option==3):
            self.valid_lambda_array = Valid_Dict['V_LAM']
            self.valid_reflectivity_array = Valid_Dict['V_REF']
            self.valid_transmissivity_array = Valid_Dict['V_TRANS']
            self.valid_emissivity_array = Valid_Dict['V_EMISS']
        elif (self.validation_option==2):
            self.valid_theta_array = Valid_Dict['V_THETA']
            self.valid_ref_vs_theta = Valid_Dict['V_REF_V_THETA']
            self.valid_trans_vs_theta = Valid_Dict['V_TRANS_V_THETA']
            self.valid_emiss_vs_theta = Valid_Dict['V_EMISS_V_THETA']
            
        return 1
            
        #self.valid_lambda_array = []
        #self.valid_reflectivity_array = []
        #self.valid_transmissivity_array = []
        #self.valid_emissivity_array = []
        #self.vald_theta_array = []
        #self.valid_ref_vs_theta = []
        #self.valid_trans_vs_theta = []
        #self.valid_emiss_vs_theta = []
        
    ### Methods to compute all Fresnel quantities at once!
    ### to compute the emissivity, one needs to compute R and T anyway
    ### so  might as well compute them all at once
    def fresnel(self):
        nc = np.zeros(len(self.d),dtype=complex)
        for i in range(0,len(self.lambda_array)):
            for j in range(0,len(self.d)):
                nc[j] = self.n[j][i]
                
            k0 = np.pi*2/self.lambda_array[i]
            ### get transfer matrix for this k0, th, pol, nc, and d
            M = tmm.tmm(k0, self.theta, self.pol, nc, self.d)
            ### get t amplitude
            t = 1./M["M11"]
            ### get incident/final angle
            ti = M["theta_i"]
            tL = M["theta_L"]
            ### get geometric factor associated with transmission
            fac = nc[len(self.d)-1]*np.cos(tL)/(nc[0]*np.cos(ti))
            ### get reflection amplitude
            r = M["M21"]/M["M11"]
            ### get Reflectivity
            self.reflectivity_array[i] = np.real(r * np.conj(r))
            ### get Transmissivity
            self.transmissivity_array[i] = np.real(t*np.conj(t)*fac)
            self.emissivity_array[i] = 1 - self.reflectivity_array[i] - self.transmissivity_array[i]

        return 1
    ### currently will return derivative of reflectivity and emissivity 
    ### wrt to thickness of layer i
    def fresnel_prime(self, layer_i):
        nc = np.zeros(len(self.d),dtype=complex)
        for i in range(0,len(self.lambda_array)):
            for j in range(0,len(self.d)):
                nc[j] = self.n[j][i]
                
            k0 = np.pi*2/self.lambda_array[i]
            ### get transfer matrix for this k0, th, pol, nc, and d
            M = tmm.tmm(k0, self.theta, self.pol, nc, self.d)
            dM_ds = tmm.d_tmm_dsi(layer_i, k0, self.theta, self.pol, nc, self.d)
            
            ### store all relevant matrix elements in variable names
            M21 = M["M21"]
            M11 = M["M11"]
            M21p = dM_ds["dM21_ds"]
            M11p = dM_ds["dM11_ds"]
            
            ### get reflection amplitude
            r = M["M21"]/M["M11"]
            r_star = np.conj(r)
            
            ### get derivative of reflection amplitudes
            r_prime = (M11*M21p - M21*M11p)/(M11*M11)
            r_prime_star = np.conj(r_prime)
            
            ### get reflectivity
            R_prime = r_prime * r_star + r * r_prime_star
            ### get Reflectivity
            self.reflectivity_prime_array[i] = np.real(R_prime)
            ### get Transmissivity
            self.emissivity_prime_array[i] = 1 - self.reflectivity_prime_array[i] 

        return 1
    
    def angular_fresnel(self, lambda_0):
        ### create an array for RI of each layer at the
        ### desired wavelength
        nc = np.zeros(len(self.d),dtype=complex)
        ### get RI for each layer and store it in nc array
        #for i in range(0,len(self.matlist)):
        #    nc[i] = datalib.Material_RI(lambda_0, self.matlist[i])
        idx, = np.where(self.lambda_array <= lambda_0)
        idx_val = idx[len(idx)-1]
        for i in range(0,len(self.matlist)):
            nc[i] = self.n[i][idx_val]
            
        k0 = np.pi*2/lambda_0
        i=0
        for thetai in self.theta_array:
            ### increment by 1/2 degrees
            M = tmm.tmm(k0, thetai, self.pol, nc, self.d)
            
            t = 1./M["M11"]
            
            ### get incident/final angle
            ti = M["theta_i"]
            tL = M["theta_L"]
            ### get geometric factor associated with transmission
            fac = nc[len(self.d)-1]*np.cos(tL)/(nc[0]*np.cos(ti))
            ### get reflection amplitude
            r = M["M21"]/M["M11"]
            ### get Reflectivity
            self.r_vs_theta[i] = np.real(r * np.conj(r))
            ### get Transmissivity
            self.t_vs_theta[i] = np.real(t*np.conj(t)*fac)
            self.eps_vs_theta[i] = 1 - self.r_vs_theta[i] - self.t_vs_theta[i]
            i = i+1
            
        return 1

    ### In case users ONLY wants reflectivity
    def reflectivity(self):
        nc = np.zeros(len(self.d),dtype=complex)
        for i in range(0,len(self.lambda_array)):
            for j in range(0,len(self.d)):
                nc[j] = self.n[j][i]
                
            k0 = np.pi*2/self.lambda_array[i]
            self.reflectivity_array[i] = tmm.Reflect(k0, self.theta, self.pol, nc, self.d)

        return 1
    ### In case user ONLY wants transmissivity
    def transmissivity(self):
        nc = np.zeros(len(self.d),dtype=complex)
        for i in range(0,len(self.lambda_array)):
            for j in range(0,len(self.d)):
                nc[j] = self.n[j][i]
                
            k0 = np.pi*2/self.lambda_array[i]
            self.transmissivity_array[i] = tmm.Trans(k0, self.theta, self.pol, nc, self.d)

        return 1
    
    ### Fresnel methods when explicit angle-averaging is requested...
    ### Need to have FOM methods to accompany this
    def fresnel_ea(self):
        if (self.explicit_angle!=1):
            error = 'ERROR: EXPLIT ANGLE OPTION NOT SELECTED! \n'
            error = error + 'RE-INSTANTIATE YOUR MULTILAYER CLASS AND BE SURE \n'
            error = error + 'TO INCLUDE A LINE IN YOUR STRUCTURE DICTIONARY LIKE THE FOLLOWING: \n'
            error = error + 'EXPLICIT_ANGLE: 1'
            print(error)
            exit()

        ### The angles come from Gauss-Legendre quadrature
        nc = np.zeros(len(self.d),dtype=complex)
        ### outter loop is over wavelength - this modulates the RI
        for i in range(0,len(self.lambda_array)):
            k0 = np.pi*2/self.lambda_array[i]
            ### for given wavelength, the stack will have the following set of RIs
            for j in range(0,len(self.d)):
                nc[j] = self.n[j][i]
            
            ### iterate over angles
            for j in range(0,len(self.t)):
                ### for given angle, k0, pol, nc, and d
                ### compute M
                Mp = tmm.tmm(k0, self.t[j], 'p', nc, self.d)
                Ms = tmm.tmm(k0, self.t[j], 's', nc, self.d)
                ### get amplitudes and related quantities from M dictionaries
                tp = 1./Mp["M11"]
                ts = 1./Ms["M11"]
                tp_i = Mp["theta_i"]
                ts_i = Ms["theta_L"]
                tp_L = Mp["theta_L"]
                ts_L = Ms["theta_L"]
                facp = nc[len(self.d)-1]*np.cos(tp_L)/(nc[0]*np.cos(tp_i))
                facs = nc[len(self.d)-1]*np.cos(ts_L)/(nc[0]*np.cos(ts_i))
                rp = Mp["M21"]/Mp["M11"]
                rs = Ms["M21"]/Ms["M11"]
                
                ### Reflectivity for each polarization
                self.reflectivity_array_p[j][i] = np.real(rp * np.conj(rp))
                self.reflectivity_array_s[j][i] = np.real(rs * np.conj(rs))
                ### Transmissivity for each polarization
                self.transmissivity_array_p[j][i] = np.real(tp*np.conj(tp)*facp)
                self.transmissivity_array_s[j][i] = np.real(ts*np.conj(ts)*facs)
                ### Emissivity for each polarization
                self.emissivity_array_p[j][i] = 1. - self.reflectivity_array_p[j][i] - self.transmissivity_array_p[j][i]
                self.emissivity_array_s[j][i] = 1. - self.reflectivity_array_s[j][i] - self.transmissivity_array_s[j][i]
                
        return 1
    
    ### Method to evaluate/update thermal emission spectrum - normal angle only!
    def thermal_emission(self):
        ### Temperature might change, update BB spectrum
        self.BBs = datalib.BB(self.lambda_array, self.T_ml)
        ### Emissivity doesn't change unless structure changes
        self.thermal_emission_array = self.BBs * self.emissivity_array
        return 1
    
    ### Method to evaluate/update thermal emission spectrum
    def thermal_emission_ea(self):
        
        ### Temperature might change, update BB spectrum
        self.BBs = datalib.BB(self.lambda_array, self.T_ml)
        #temp = np.zeros(len(self.lambda_array))
        
        for i in range(0,len(self.t)):
            ### Thermal emission goes like BBs(lambda) * eps(theta, lambda) * cos(theta)
            for j in range(0,len(self.lambda_array)):
                self.thermal_emission_array_p[i][j] = self.BBs[j] * self.emissivity_array_p[i][j] * np.cos(self.t[i])
                self.thermal_emission_array_s[i][j] = self.BBs[j] * self.emissivity_array_s[i][j] * np.cos(self.t[i])
            
        return 1
    
    ''' METHOD FOR J-AGG ENHANCEMENT '''
    def jagg_sd(self):
        self.jagg_sd_val = numlib.Integrate(self.emissivity_array, self.lambda_array, 500e-9, 700e-9)/200e-9
        return 1
    
    ''' METHODS FOR STPVLIB!!! '''
    
    ### Normal versions first - no explicit dependence on angle
    
    ### Spectral Efficiency - see Eq. 4 in Jeon et al, Adv. Energy Mater. 2018 (8) 1801035
    def stpv_se(self):
        self.spectral_efficiency_val = stpvlib.SpectralEfficiency(self.thermal_emission_array, self.lambda_array, self.lbg)
        return 1
    
    ### Power density - see Eq. 3 in Jeon et al, Adv. Energy Mater. 2018 (8) 1801035
    def stpv_pd(self):
        self.power_density_val = stpvlib.Pwr_den(self.thermal_emission_array, self.lambda_array, self.lbg)
        return 1
    
    ### TPV Efficiency, see Eq. S20-S26 in Jeon et al, Adv. Energy Mater. 2018 (8) 1801035
    def stpv_etatpv(self):
        self.tpv_efficiency_val = stpvlib.Eta_TPV(self.thermal_emission_array, self.lambda_array, self.PV, self.T_cell)
        return 1
    
    ### Explicit Angle versions of methods for STPV quantities
    def stpv_se_ea(self):
        self.spectral_efficiency_val = stpvlib.SpectralEfficiency_EA(self.thermal_emission_array_p, self.thermal_emission_array_s, self.lambda_array, self.lbg, self.t, self.w)

        
    def stpv_pd_ea(self):
        self.power_density_val = stpvlib.Pwr_den_EA(self.thermal_emission_array_p, self.thermal_emission_array_s, self.lambda_array, self.lbg, self.t, self.w)
        return 1
    
    def stpv_etatpv_ea(self):
        self.tpv_efficiency_val = stpvlib.Eta_TPV_EA(self.thermal_emission_array_p, self.thermal_emission_array_s, self.lambda_array, self.PV, self.T_cell, self.t, self.w)
    
    
    ### Absorber Efficiency - see 
    def stpv_etaabs(self):
        alpha = stpvlib.absorbed_power_ea(self.lambda_array, self.n, self.d, self.solarconc)
        beta = stpvlib.p_in(self.thermal_emission_array, self.lambda_array)
        print("alpha is ",alpha)
        print("beta is ",beta)
        self.absorber_efficiency_val = (alpha - beta)/alpha
        return 1
        
    def stpv_etaabs_ea(self):
        
        ### Power absorbed is going to explicitly consider a range of incident angles which
        ### will depend on the solar concentration
        alpha = stpvlib.absorbed_power_ea(self.lambda_array, self.n, self.d, self.solarconc)
        beta = stpvlib.p_in_ea(self.thermal_emission_array_p, self.thermal_emission_array_s, self.lambda_array, self.t, self.w )
        self.absorber_efficiency_val = (alpha - beta)/alpha
        return 1
    
    ### Currently using optimistic values for Voc (0.706 mV) and FF (0.828) reported
    ### in Nanoscale Research Letters, (2016) vol 11 pg 453, L.-X. Wang, Z.-Q. Zhou, T.-N. Zhang, X. Chen
    ### and M. Lu
    def pv_conversion_efficiency(self):
        self.short_circuit_current_val = stpvlib.ambient_jsc(self.emissivity_array, self.lambda_array, self.lbg)
        #self.open_circuit_voltage_val = stpvlib.Voc(self.short_circuit_current_val, self.T_cell)
        #self.fill_factor_val = stpvlib.FF(self.open_circuit_voltage_val, self.T_cell)
        self.incident_power =  stpvlib.integrated_solar_power(self.lambda_array)
        self.conversion_efficiency_val = self.short_circuit_current_val * 0.828 * 0.706
        #self.conversion_efficiency_val = self.conversion_efficiency_val*self.open_circuit_voltage_val
        #self.conversion_efficiency_val = self.conversion_efficiency_val*self.fill_factor_val
        self.conversion_efficiency_val = self.conversion_efficiency_val / self.incident_power

        return 1        

    def step_emissivity(self, lambda_0, delta_lambda):
        idx = 0
        for lam in self.lambda_array:
            if (lam > (lambda_0 - delta_lambda/2) and lam < (lambda_0 + delta_lambda/2)):
                self.emissivity_array[idx] = 1.
                self.transmissivity_array[idx] = 0.
                self.reflectivity_array[idx] = 0.
                idx = idx + 1
            else:
                self.emissivity_array[idx] = 0.
                self.transmissivity_array[idx] = 0.
                self.reflectivity_array[idx] = 1.
                idx = idx + 1
        return 1

    def step_emissivity_ea(self, lambda_0, delta_lambda):
        ### The angles come from Gauss-Legendre quadrature
        ### outter loop is over wavelength - this modulates the RI
        for j in range(0, len(self.t)):
            idx = -1
            for lam in self.lambda_array:
                idx = idx + 1
                if (lam>(lambda_0 - delta_lambda/2) and lam < (lambda_0 + delta_lambda/2)):
                    self.emissivity_array_p[j][idx] = 1.
                    self.emissivity_array_s[j][idx] = 1.
                    self.reflectivity_array_p[j][idx] = 0.
                    self.reflectivity_array_s[j][idx] = 0.
                    self.transmissivity_array_p[j][idx] = 0.
                    self.transmissivity_array_s[j][idx] = 0.
                else:
                    self.emissivity_array_p[j][idx] = 0.
                    self.emissivity_array_s[j][idx] = 0.
                    self.reflectivity_array_p[j][idx] = 1.
                    self.reflectivity_array_s[j][idx] = 1.
                    self.transmissivity_array_p[j][idx] = 0.
                    self.transmissivity_array_s[j][idx] = 0.
        return 1
                    

    def step_reflectivity(self, lambda_0, delta_lambda):
        idx = 0
        for lam in self.lambda_array:
            if (lam > (lambda_0 - delta_lambda/2) and lam < (lambda_0 + delta_lambda/2)):
                self.emissivity_array[idx] = 0.
                self.transmissivity_array[idx] = 0.
                self.reflectivity_array[idx] = 1.
                idx = idx + 1
            else:
                self.emissivity_array[idx] = 1.
                self.transmissivity_array[idx] = 0.
                self.reflectivity_array[idx] = 0.
                idx = idx + 1
        return 1
    
                
    ''' METHODS FOR COLORLIB!!! '''
    
    ### displays the percieved color of an object at a specific temperature
    ### based only on thermal emission
    def thermal_color(self):
        string = "Color at T = " + str(self.T_ml) + " K"
        colorlib.RenderColor(self.thermal_emission_array, self.lambda_array, string)
        self.thermal_rgb = colorlib.RGB_FromSpec(self.thermal_emission_array, self.lambda_array)
        return 1
    
    ### Displays the perceived color of an object based only
    ### on reflected light
    def ambient_color(self):
        string = "Ambient Color"
        colorlib.RenderColor(self.reflectivity_array, self.lambda_array, string)
        self.reflective_rgb = colorlib.RGB_FromSpec(self.reflectivity_array, self.lambda_array)
        return 1
    
    ### Displays the percieved color of a narrow bandwidth lightsource
    def pure_color(self, wl):
        Spectrum = np.zeros_like(self.lambda_array)
        for i in range(0,len(Spectrum)):
            if abs(self.lambda_array[i] - wl)<5e-9:
                Spectrum[i] = 1
        colorlib.RenderColor(Spectrum, self.lambda_array, str(wl))
        return 1
    
    def classify_color(self):
        self.color_name = colorlib.classify_color(self.reflectivity_array, self.lambda_array)
    
    ''' METHODS FOR LIGHTLIB '''
    def luminous_efficiency(self):
        self.luminous_efficiency_val = lightlib.Lum_efficiency(self.lambda_array, self.thermal_emission_array)
        return 1
    def luminous_efficiency_prime(self):
        self.luminous_efficiency_prime_val = lightlib.Lum_efficiency_prime(self.lambda_array, self.thermal_emission_array, self.BBs*self.emissivity_prime_array)
    
    def normalized_luminous_power(self):
        self.luminous_power_val = lightlib.normalized_power(self.lambda_array, self.thermal_emission_array, self.BBs)
        return 1
    
    ''' METHODS FOR COOLINGLIB !!! '''
    def cooling_power(self):
        self.radiative_power_val = coolinglib.Prad(self.thermal_emission_array_p, self.thermal_emission_array_s, self.lambda_array, self.t, self.w)
        self.atmospheric_power_val = coolinglib.Patm(self.emissivity_array_p, self.emissivity_array_s, self.T_amb, self.lambda_array, self.t, self.w)
        self.solar_power_val = coolinglib.Psun(self.theta_sun, self.lambda_array, self.n, self.d)
        self.cooling_power_val = self.radiative_power_val - self.atmospheric_power_val - self.solar_power_val
        return 1
    
    ''' MISCELLANEOUS METHODS TO MANIPULATE THE STRUCTURE
        OR GATHER DATA ABOUT THE STRUCTURE '''
    
    ### Method to add a layer to the bottom of the structure
    ### and re-compute desired quantities
    def insert_layer(self, layer_number, material, thickness):
        ### just use numpy insert for thickness array
        new_d = np.insert(self.d, layer_number, thickness)
        ### because material names can have variable number of characters,
        ### insert may not reliably work... do "manually"
        new_m = []
        for i in range(0,layer_number):
            new_m.append(self.matlist[i])
        new_m.append(material)
        for i in range(layer_number+1,len(self.matlist)+1):
            new_m.append(self.matlist[i-1])

        ### de-allocate memory associated with self.d, self.matlist, self.n arrays
        self.d = None
        self.matlist = None
        self.n = None

        ### assign new values to self.d, self.matlist, self.n
        self.d = new_d
        self.matlist = new_m
  
        self.n = None 
        self.n = np.zeros((len(self.d),len(self.lambda_array)),dtype=complex)
        for i in range(0,len(self.matlist)):
                self.n[:][i] = datalib.Material_RI(self.lambda_array, self.matlist[i])
        
        ### in all cases, updated Fresnel quantities
        self.fresnel()
        
        if (self.explicit_angle):
            
            self.fresnel_ea()
        
        ### if a thermal application is requested, update Thermal Emission as well
        #if (self.stpv_emitter_calc or self.stpv_absorber_calc or self.cooling_calc or self.lightbulb_calc or self.color_calc):
        #    self.ThermalEmission()
        if self.stpv_emitter_calc:
            self.thermal_emission()
            self.stpv_se()
            self.stpv_pd()
            self.stpv_etatpv()
            
            if (self.explicit_angle):
                self.thermal_emission_ea()
                self.stpv_se_ea()
                self.stpv_pd_ea()
                ### need to implement eta_tpv_ea method.
                #self.stpv_etatpv_ea()
            
        if self.stpv_absorber_calc:
            self.thermal_emission()
            
            if (self.explicit_angle):
                self.thermal_emission_ea()
                self.stpv_etaabs()
            else:
                self.stpv_etaabs()
        
        
        return 1

    ### Get the RI of a particular layer (at each wavelength specified by user)
    def layer_ri(self, layer):
        RI = np.zeros(len(self.lambda_array),dtype=complex)
        for i in range(0,len(self.lambda_array)):
            RI[i] = self.n[layer][i]
        return RI
    
    ### Define the RI of a specified layer to be an alloy
    ### between two specified materials, mat1 and mat2,
    ### using Bruggenmans approximation
    def layer_alloy(self, layer, fraction, mat1, mat2, model):
        ### Bruggeman model must be specified
        if (model=='Bruggeman'):
            ### Get RIs of two materials... mat1 can be 
            ### a string that codes a material name or 
            ### it can be a single number
            if(isinstance(mat1, str)):
                n_1 = datalib.Material_RI(self.lambda_array, mat1)
            else:
                n_1 = mat1
                
            n_2 = datalib.Material_RI(self.lambda_array, mat2)
            
            for i in range(0,len(self.lambda_array)):
                if(isinstance(mat1, str)):
                    eps1 = n_1[i]*n_1[i]
                else:
                    eps1 = n_1*n_1
                    
                eps2 = n_2[i]*n_2[i]
                flag = 1
                f1 = (1-fraction)
                f2 = fraction
                b = (2*f1-f2)*eps1 + (2*f2 - f1)*eps2
                arg = 8*eps1*eps2 + b*b
                srarg = np.sqrt(arg)
                
                if (np.imag(arg)<0):
                    flag = -1
                else:
                    flag = 1
                    
                epsBG = (b+flag*srarg)/4.
                self.n[layer][i] = np.sqrt(epsBG)
        #### Default is Maxwell-Garnett        
        else:
            if(isinstance(mat1, str)):
                n_1 = datalib.Material_RI(self.lambda_array, mat1)
            else:
                n_1 = mat1
                
            n_2 = datalib.Material_RI(self.lambda_array, mat2)
            f = fraction
            

            for i in range(0,len(self.lambda_array)):
                ### eps1 == epsD and eps2 == epsM in MG notation
                if(isinstance(mat1, str)):
                    epsD = n_1[i]*n_1[i]
                else:
                    epsD = n_1*n_1

                epsM = n_2[i]*n_2[i]
                num = epsD*(2*f*(epsM-epsD) + epsM + 2*epsD)
                denom = 2*epsD + epsM + f*(epsD-epsM)
                self.n[layer][i] = np.sqrt((num/denom))
                
        return 1
    
    #### Sets the refractive index for a specified layer
    #### to a single specified refractive index value
    def layer_static_ri(self, layer, RI):
        for i in range(0,len(self.lambda_array)):
            self.n[layer][i] = RI
            
        return 1

    def layer_lorentz(self, layer, omega_p, omega_0, gamma):
        c = 299792458.
        ci = 0+1j
        for i in range(0,len(self.lambda_array)):
            omega = 2*np.pi*c/self.lambda_array[i]
            eps_lr = 1 + omega_p**2/(omega_0**2 - omega**2 - ci*omega*gamma)
            self.n[layer][i] = np.sqrt(eps_lr)
        return 1
    ### METHODS FOR PLOTTING DATA!
    
    ### Plot thermal emission
    def plot_te(self):
        plt.plot(self.lambda_array*1e9, self.thermal_emission_array, 'red')
        string = "Thermal Emission at " + str(self.T_ml) + " K"
        plt.legend(string)
        plt.show()
        return 1
    
    ### Plot reflectivity
    def plot_reflectivity(self):
        plt.plot(self.lambda_array*1e9, self.reflectivity_array, 'red')
        string = "Reflectivity"
        plt.legend(string)
        plt.show()
        return 1
    
    ### Plot emissivity
    def plot_emissivity(self):
        plt.plot(self.lambda_array*1e9, self.emissivity_array, 'blue')
        string = "Emissivity"
        plt.legend(("Emissivity"))
        plt.show()
        return 1   
    
    ### RESONANCE METHODS
    
    ### Find SPP mode for a structure at a particular
    ### wavelength that has index idx in the array of wavelengths!

    def find_spp(self, idx):
        ### get wavevector at idx'th wavelength
        k0 = np.pi*2/self.lambda_array[idx]
        L = len(self.d)
        ### array of RIs for the structure at the idx'th wavelength
        nc = np.zeros(L,dtype=complex)
        for j in range(0,L):
            nc[j] = self.n[j][idx]
            
        ### SPP is "above the light line" so
        ### start searching beta at the light line
        b_beg = k0*nc[L-1]
        ### set maximum beta as k0*index of incident material
        b_end = k0*nc[0]
        
        ### alpha is harder to bound... must be positive, should be much smaller
        ### than beta
        a_beg = 0.000001
        a_end = 0.2*b_end
        
        beta = np.linspace(b_beg, b_end,100)
        alpha = np.linspace(a_beg, a_end, 100)
        
        ### initialize values
        rr_max = -100
        rr_temp = 0
        a_spp = 0
        b_spp = 0

        for a in alpha:
           # print(" ")
            for b in beta:
                kx = b + a*1j
                rr_temp = tmm.tmm_ab(k0, kx, 'p', nc, self.d)                
                #print(np.real(kx),np.imag(kx), t_temp)
                #t_array[k] = t_temp
                #k+=1
                if rr_temp>rr_max:
                    rr_max = rr_temp
                    a_spp = a
                    b_spp = b

        self.spp_resonance_val = b_spp+a_spp*1j
        return 1

    def find_pa(self, idx):
        ### get wavevector at idx'th wavelength
        k0 = np.pi*2/self.lambda_array[idx]
        L = len(self.d)
        ### array of RIs for the structure at the idx'th wavelength
        nc = np.zeros(L,dtype=complex)
        for j in range(0,L):
            nc[j] = self.n[j][idx]
            
        ### SPP is "above the light line" so
        ### start searching beta at the light line
        b_beg = k0*nc[L-1]
        ### set maximum beta as k0*index of incident material
        b_end = k0*nc[0]
        
        ### alpha is harder to bound... must be positive, should be much smaller
        ### than beta
        a_beg = 0.000001
        a_end = 0.2*b_end
        
        beta = np.linspace(b_beg, b_end,100)
        alpha = np.linspace(a_beg, a_end, 100)
        
        ### initialize values
        rr_min = 100
        rr_temp = 0
        a_spp = 0
        b_spp = 0

        for a in alpha:
           # print(" ")
            for b in beta:
                kx = b + a*1j
                rr_temp = tmm.tmm_ab(k0, kx, 'p', nc, self.d)                
                #print(np.real(kx),np.imag(kx), t_temp)
                #t_array[k] = t_temp
                #k+=1
                if rr_temp<rr_min:
                    rr_min = rr_temp
                    a_spp = a
                    b_spp = b

        self.pa_resonance_val = b_spp+a_spp*1j
        return 1

    
    def inline_structure(self, args):
        if 'Lambda_List' in args:
            lamlist = args['Lambda_List']
            self.lambda_array = np.linspace(lamlist[0],lamlist[1],int(lamlist[2]))
        else:
            print(" Lambda array not specified! ")
            print(" Choosing default array of 1000 wl between 400 and 6000 nm")
            self.lambda_array = np.linspace(400e-9,6000e-9,1000)


        
        if 'Thickness_List' in args:
            self.d = args['Thickness_List']
        ### default structure
        else:
            print("  Thickness array not specified!")
            print("  Proceeding with default structure - optically thick W! ")
            self.d = [0, 900e-9, 0]
            self.matlist = ['Air', 'W', 'Air']
            self.n = np.zeros((len(self.d),len(self.lambda_array)),dtype=complex)
            for i in range(0,len(self.matlist)):
                self.n[:][i] = datalib.Material_RI(self.lambda_array, self.matlist[i])
                
        if 'Material_List' in args:
            self.matlist = args['Material_List']
            self.n = np.zeros((len(self.d),len(self.lambda_array)),dtype=complex)
            for i in range(0,len(self.matlist)):
                    self.n[:][i] = datalib.Material_RI(self.lambda_array, self.matlist[i])
            
        else:
            print("  Material array not specified!")
            print("  Proceeding with default structure - optically thick W! ")
            self.d = [0, 900e-9, 0]
            self.matlist = ['Air', 'W', 'Air']
            self.n = np.zeros((len(self.d),len(self.lambda_array)),dtype=complex)
            for i in range(0,len(self.matlist)):
                    self.n[:][i] = datalib.Material_RI(self.lambda_array, self.matlist[i])

        ### Temperature arguments!
        if 'Temperature' in args:
            self.T_ml = args['Temperature']
        elif 'Structure_Temperature' in args:
            self.T_ml = args['Structure_Temperature']
        else:
            print(" Temperature not specified!")
            print(" Proceeding with default T = 300 K")
            self.T_ml = 300
        if 'PV_Temperature' in args:
            self.T_cell = args['PV_Temperature']
        else:
            self.T_cell = 300
        if 'Ambient_Temperature' in args:
            self.T_amb = args['Ambient_Temperature']
        else:
            self.T_amb = 300
        
        ### Check to see what calculations should be done!
        if 'STPV_EMIT' in args:
            self.stpv_emitter_calc = args['STPV_EMIT']
        else:
            self.stpv_emitter_calc = 0
        if 'STPV_ABS' in args:
            self.stpv_absorber_calc = args['STPV_ABS']
        else:
            self.stpv_absorber_calc = 0
        if 'COOLING' in args:
            self.cooling_calc = args['COOLING']
        else:
            self.cooling_calc = 0
        if 'LIGHTBULB' in args:
            self.lightbulb_calc = args['LIGHTBULB']
        else:
            self.lightbulb_calc = 0
        if 'COLOR' in args:
            self.color_calc = args['COLOR']
        else:
            self.color_calc = 0
        if 'EXPLICIT_ANGLE' in args:
            self.explicit_angle = args['EXPLICIT_ANGLE']
        else:
            self.explicit_angle = 0
        if 'DEG' in args:
            self.deg = args['DEG']
        else:
            self.deg = 7
            
        return 1

    
