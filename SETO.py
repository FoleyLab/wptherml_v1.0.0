### Import WPTHERML class!
from wptherml.wpml import multilayer
from matplotlib import pyplot as plt
from wptherml.datalib import datalib
import numpy as np

### Define structure!
### Structure II in Fig 2
structureII = {

        ### will modify layer 2 to be an alloy later, for now it is HfO2
        'Material_List': ['Air', 'Al2O3', 'HfO2', 'W', 'Air'],
        'Thickness_List': [0, 63e-9, 94e-9, 900e-9, 0],
        'Lambda_List': [400e-9, 6000e-9, 2000],
        'Temperature': 1023,
        'EXPLICIT_ANGLE': 0,
        'STPV_EMIT': 1
     
        }
### create instance of structureII called sII
sII = multilayer(structureII)

### change layer 2 of sII to be the 17% TiN in HfO2 alloy
sII.layer_alloy(2, 0.17, 'HfO2', 'TiN', 'Bruggeman')

### need to recompute otical properties and figures of merit with the alloy layer
sII.fresnel()
#sII.fresnel_ea()

py = datalib.Abs_Pyromark(sII.lambda_array)
BBs = datalib.BB(sII.lambda_array, sII.T_ml)
AMs = datalib.AM(sII.lambda_array)

plt.plot(sII.lambda_array*1e9, 600*sII.emissivity_array*AMs/1e9/1e3, sII.lambda_array*1e9, 600*py*AMs/1e9/1e3)
plt.plot(sII.lambda_array*1e9, 2*np.pi*sII.thermal_emission_array/1e9/1e3, sII.lambda_array*1e9, 2*np.pi*py*BBs/1e9/1e3, sII.lambda_array*1e9, 2*np.pi*BBs/1e9/1e3)
plt.show()


'''
T_List = np.linspace(800, 1400, 61)

for Tmp in T_List:
    sII.T_ml = Tmp+273
    sII.thermal_emission()
    sII.thermal_emission_ea()
    ### recompute figures of merit, using explicit angle dependence of emissivity
    sII.stpv_se_ea()
    sII.stpv_pd_ea()
    sII.stpv_etatpv_ea()
    print(sII.T_ml,sII.spectral_efficiency_val,1e-4*sII.power_density_val, sII.tpv_efficiency_val)

'''