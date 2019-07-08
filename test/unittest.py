from wptherml.wpml import multilayer
import numpy as np

''' The following ml structure uses all materials supported as of 06/30/2019
    and systematically tests the geometry and material parsing, computation of optical properties,
    computational of thermal emission properties, and computation of figures of merit '''
structure = {
  'Material_List': ['Air', 'HfO2', 'Al2O3', 'SiO2', 'AlN', 'W', 'W', 'Re', 'Rh', 'Ru', 'Ag', 'Au', 'Pd', 'Pt', 'Si', 'TiO2',  'Air'],
  'Thickness_List': [0, 20e-9, 21e-9, 22e-9, 23e-9, 24e-9, 25e-9, 26e-9, 27e-9, 28e-9, 29e-9, 30e-9, 31e-9, 32e-9, 33e-9, 34e-9,  0],
  'Lambda_List': [500e-9, 600e-9, 101]

}
ml = multilayer(structure)

''' test to make sure thicknesses of each layer are parsed correctly'''
def test_thickness_parsing():
    assert( ml.d[1] == 20e-9)
    assert( ml.d[2] == 21e-9)
    assert( ml.d[3] == 22e-9)
    assert( ml.d[4] == 23e-9)
    assert( ml.d[5] == 24e-9)
    assert( ml.d[6] == 25e-9)
    assert( ml.d[7] == 26e-9)
    assert( ml.d[8] == 27e-9)
    assert( ml.d[9] == 28e-9)
    assert( ml.d[10] == 29e-9)
    assert( ml.d[11] == 30e-9)
    assert( ml.d[12] == 31e-9)
    assert( ml.d[13] == 32e-9)
    assert( ml.d[14] == 33e-9)
    assert( ml.d[15] == 34e-9)

''' test to make sure RI of each material parsed correctly, should
    compare to RI of each material in the folloing list at 501 nm:
'Air', 'HfO2', 'Al2O3', 'SiO2', 'AlN', 'W', 'W', 'Re', 'Rh', 'Ru', 'Ag', 'Au', 'Pd', 'Pt', 'Si', 'TiO2',  'Air' '''
def test_ri_parsing():
    ri_array = np.array([ 1.00000000 +0.00000000e+00j,  1.97985108 +0.00000000e+00j,
  1.73984208 +0.00000000e+00j,  1.48985208 +0.00000000e+00j,
  3.21397468 +0.00000000e+00j,  3.40033724 +2.69524924e+00j,
  3.40033724 +2.69524924e+00j,  4.01405239 +2.64618165e+00j,
  1.88512950 +4.68073102e+00j,  3.30699692 +4.79672201e+00j,
  0.04997774 +3.14242568e+00j,  0.94744365 +1.87688826e+00j,
  1.53073694 +3.57015516e+00j,  1.97802044 +3.45088214e+00j,
  4.28360000 +4.80709000e-02j,  2.47960735 +6.15938472e-07j,
  1.00000000 +0.00000000e+00j],dtype=complex)
    
    for i in range(0,len(ri_array)):
        assert( np.abs( ml.n[i,1] -  ri_array[i] ) < 1e-6 )

''' test to make sure normal-incidence fresnel calculations are working properly... 
    Will test the reflectivity, transmissivity, and emissivity of the above multilayer stack at 501 nm '''
def test_normal_fresnel():
    ### test values at 501 nm
    assert( np.abs(ml.reflectivity_array[1] - 0.126899985878) < 1e-6)
    assert( np.abs(ml.emissivity_array[1] - 0.873100013379 ) < 1e-6)
    assert( np.abs(ml.transmissivity_array[1] - 7.43394762387e-10) < 1e-6)
    ### test sum of values at 510 nm... should be 1 within round off error
    total = ml.reflectivity_array[10] + ml.emissivity_array[10] + ml.transmissivity_array[10]
    assert( np.abs( total - 1.0) < 1e-6)
    
''' test to make sure normal thermal emission calculation is orking properly...
    will test the normal thermal emission at T = 1000 K and lambda = 501 nm '''
def test_normal_thermal_emission():
    ml.T_ml = 1000
    ml.thermal_emission()
    assert( np.abs(ml.thermal_emission_array[1] - 1110.92749199) < 1e-6)
  
''' test to make sure normal STPV methods are working propertly...
    will test the the spectral efficiency, power density, and tpv efficiency
    for the above multilayer at 1000 K '''
def test_normal_stpvlib():
    ml.stpv_se()
    ml.stpv_pd()
    ml.stpv_etatpv()
    assert( np.abs(ml.spectral_efficiency_val -  0.255520461044) < 1e-6)
    assert( np.abs(ml.power_density_val - 0.00119731168356) < 1e-6)
    assert( np.abs(ml.tpv_efficiency_val - 0.00185689588813) < 1e-6)

''' test to make sure normal LIGHTLIB methods are working propertly...
    will test the luminous efficiency of the multilayer at 1000 K '''
def test_normal_lightlib():
    ml.luminous_efficiency()
    assert( np.abs(ml.luminous_efficiency_val - 0.863383797679) < 1e-6)

structure['EXPLICIT_ANGLE'] = 1
ml = multilayer(structure)
    
''' test to make sure explicit-angle fresnel methods are working properly...
    will test the p- and s-polarized reflectivity, transmissivity, and emissivity 
    at 501 nm and 45 deg '''
def test_explicit_angle_fresnel():
    ml.pol = 'p'
    ml.fresnel_ea()
    assert( np.abs(ml.t[3]*180/np.pi - 45.0) < 1e-6)
    assert( np.abs(ml.reflectivity_array_p[3,1] - 0.0630899230711) < 1e-6)
    assert( np.abs(ml.emissivity_array_p[3,1] - 0.936910076294) < 1e-6)
    assert( np.abs(ml.transmissivity_array_p[3,1] - 6.35171792001e-10) < 1e-6)
    assert( np.abs(ml.reflectivity_array_s[3,1] - 0.166381585217) < 1e-6)
    assert( np.abs(ml.emissivity_array_s[3,1] - 0.833618414388) < 1e-6)
    assert( np.abs(ml.transmissivity_array_s[3,1] - 3.94488226834e-10) < 1e-6)

''' test to make sure explicit angle thermal emission is working properly...
    will test the p- and s-polarized TE at 501 nm and 45 deg at 1000 K.  '''
def test_explicit_angle_thermal_emission():
    ml.thermal_emission_ea()
    assert( np.abs(ml.thermal_emission_array_p[3,1] - 842.955466484) < 1e-6)
    assert( np.abs(ml.thermal_emission_array_s[3,1] - 750.022032157) < 1e-6)
   
''' test to make sure explicit angle STPV methods are working...
    will test e-a calculations of multilayer at 1000 K '''
def test_explicit_angle_stpvlib():
    ml.stpv_se_ea()
    ml.stpv_pd_ea()
    ml.stpv_etatpv_ea()
    ml.stpv_etaabs_ea()
    assert( np.abs(ml.spectral_efficiency_val - 0.255543236604) < 1e-6)
    assert( np.abs(ml.power_density_val - 0.00118224745587) < 1e-6)
    assert( np.abs(ml.tpv_efficiency_val - 0.00183764368281) < 1e-6)
    assert( np.abs(ml.absorber_efficiency_val - 0.999998420168) < 1e-6)
    
''' test to make sure cooling lib methods are working...
    will test multilayer at 1000 K '''    
def test_explicit_angle_coolinglib():

    ml.cooling_power()
    assert( np.abs(ml.radiative_power_val - 0.00462640871104) < 1e-6)
    assert( np.abs(ml.atmospheric_power_val - 6.80594278066e-28) < 1e-6)
    assert( np.abs(ml.solar_power_val - 120.842780408) < 1e-6)

''' test to make sure angular fresnel calculation is working properly...
    will test the reflectivity, transmissivity, and emissivity of the above multilayer stack at 
    501 nm and theta_i = 50 degrees at p- and s-polarizations '''
def test_angular_fresnel():
    ### polarization not specified, p-polarization is default
    ml.angular_fresnel(501e-9)
    assert( np.abs(ml.theta_array[100]*180/np.pi - 50) < 1e-6)
    assert( np.abs(ml.r_vs_theta[100] -  0.0510781599902) < 1e-6)
    assert( np.abs(ml.eps_vs_theta[100] - 0.948921839382) < 1e-6)
    assert( np.abs(ml.t_vs_theta[100] - 6.27663477264e-10) < 1e-6)
    ml.pol = 's'
    ml.angular_fresnel(501e-9)
    assert( np.abs(ml.r_vs_theta[100] -  0.181614746651) < 1e-6)
    assert( np.abs(ml.eps_vs_theta[100] - 0.81838525301) < 1e-6)
    assert( np.abs(ml.t_vs_theta[100] - 3.38868081285e-10) < 1e-6)
    





