

### Import WPTHERML class!
from wptherml.wpml import multilayer
from matplotlib import pyplot as plt
from wptherml.datalib import datalib
from wptherml.numlib import numlib
import numpy as np

### Define structure!
### Structure II in Fig 2
d1 = 100e-9/(2*1.76)
structureII = {

        ### will modify layer 2 to be an alloy later, for now it is HfO2
        'Material_List': ['Air', 'Al2O3', 'HfO2', 'W', 'Air'],
        'Thickness_List': [0, 60e-9, 100e-9, 900e-9, 0],
        'Lambda_List': [200e-9, 10000e-9, 5000],
        'Temperature': 1023,
        'EXPLICIT_ANGLE': 0,
        'STPV_ABS': 0
     
        }

w_only = {
        'Material_List': ['Air', 'W', 'Air'],
        'Thickness_List': [0, 900e-9, 0],
        'Lambda_List': [200e-9, 10000e-9, 5000],
        'Temperature': 1023,
        'EXPLICIT_ANGLE': 0,
        'STPV_ABS': 0
        }
### create instance of structureII called sII
sII = multilayer(structureII)
w_slab = multilayer(w_only)
w_slab.thermal_emission()

### change layer 2 of sII to be the 17% TiN in HfO2 alloy
sII.layer_alloy(2, 0.17, 'HfO2', 'TiN', 'Bruggeman')

### need to recompute otical properties and figures of merit with the alloy layer
sII.fresnel()
sII.thermal_emission()
#sII.stpv_etaabs()
#print(sII.absorber_efficiency_val)
#sII.fresnel_ea()

### pyromark absorptivity/emissivity
py = datalib.Abs_Pyromark(sII.lambda_array)
### 750 C BB spectrum
BBs = datalib.BB(sII.lambda_array, sII.T_ml)
### 600 suns
AMs = 600*datalib.AM(sII.lambda_array)
### upper bound for numerical integrals
upper = np.amax(sII.lambda_array)
### total incident solar power
solar_int = numlib.Integrate(AMs, sII.lambda_array, 1e-9, upper)
### total emitted BB power from 750 strcture
bb_int = numlib.Integrate(np.pi*BBs, sII.lambda_array, 1e-9, upper)

### pyro absorbed power
pyro_abs = numlib.Integrate(py*AMs, sII.lambda_array, 1e-9, upper)
### pyro emitted power
pyro_emit = numlib.Integrate(np.pi*BBs*py, sII.lambda_array, 1e-9, upper)

### cavity-enhanced absorbed power
ce_abs = numlib.Integrate(sII.emissivity_array*AMs, sII.lambda_array, 1e-9, upper)
### cavity-enhanced emitted power
ce_emit = numlib.Integrate(np.pi*sII.thermal_emission_array, sII.lambda_array, 1e-9, upper)

### reference w slab absorbed power
w_abs = numlib.Integrate(w_slab.emissivity_array*AMs, w_slab.lambda_array, 1e-9, upper)
### reference w slab emitted poer
w_emit = numlib.Integrate(np.pi*w_slab.thermal_emission_array, w_slab.lambda_array, 1e-9, upper)


print("  Ideal BB ")
print("Solar Power ",solar_int)
print("BB Power", bb_int)
print("Efficiency is ",100*(solar_int - bb_int)/solar_int, "%")

print("  Pyromark ")
print(" Absorbed Solar Power ", pyro_abs)
print(" Emitted Power ", pyro_emit)
print(" Efficiency ", (pyro_abs - pyro_emit)/solar_int, "%")
print(" % Absorbed ", 100*pyro_abs/solar_int)
print(" % Emitted ", 100*pyro_emit/bb_int)

print("  Cavity Enhanced ")
print(" Absorbed Solar Power ", ce_abs)
print(" Emitted Power ", ce_emit)
print(" Efficiency ", (ce_abs - ce_emit)/solar_int, "%")
print(" % Absorbed ",100*ce_abs/solar_int)
print(" % Emitted ", 100*ce_emit/bb_int)

print(" Bare W ")
print(" Absorbed Solar Power ",w_abs)
print(" Emitted Power ", w_emit)
print(" Efficiency is ",100*(w_abs - w_emit)/solar_int, "%")
print(" % Absorbed is ", 100*w_abs/solar_int)
print(" % Emitted ", 100*w_emit/bb_int)

print("Wavelength (nm), AM1.5, BB@750C, CE Emissivity, Pyromark Emissivity, W Emissivity")
for i in range(0,len(sII.lambda_array)):
  print(sII.lambda_array[i]*1e9, AMs[i]/1e12, BBs[i]/1e12, sII.emissivity_array[i], py[i], w_slab.emissivity_array[i]) 

'''
plt.plot(sII.lambda_array*1e9, AMs/1e9/1e3, 'black')
plt.plot(sII.lambda_array*1e9, py*AMs/1e9/1e3, 'red')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Power Density (W/m^2/nm)")
plt.legend(('AM 1.5 Spectrum', 'Absorbed Power'),
           loc='upper center', shadow=False)

plt.plot(sII.lambda_array*1e9, np.pi*BBs/1e9/1e3, 'black')
plt.plot(sII.lambda_array*1e9, np.pi*py*BBs/1e9/1e3, 'red')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Power Density (W/m^2/nm)")
plt.legend(('Blackbody Spectrum', 'Thermal Emission'),
           loc='upper center', shadow=False)
'''
#plt.plot(sII.lambda_array*1e9, 600*py*AMs/1e9/1e3)
#plt.show()
#plt.savefig('/Users/jay/Career/Proposals/SETO_V2/Pyro_Abs.png')


#plt.plot(sII.lambda_array*1e9, np.pi*sII.thermal_emission_array/1e9/1e3, sII.lambda_array*1e9, np.pi*py*BBs/1e9/1e3, sII.lambda_array*1e9, np.pi*BBs/1e9/1e3)
#plt.show()

