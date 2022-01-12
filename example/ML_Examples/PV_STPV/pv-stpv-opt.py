from wptherml.wpml import multilayer
from matplotlib import pyplot as plt
from wptherml.datalib import datalib
from wptherml.numlib import numlib
from scipy.optimize import minimize
import time
import numpy as np
h = 6.626e-34
c = 299792458.
temp = 428
wl1 = 3000e-9
wl2 = 4200e-9

### define Lamba_List bounds
startwl = 300e-9
endwl = 20300e-9
numwl = 5001;

### create a dummy stack to get refractive indices
structure = {
        'Material_List' : ['Air', 'SiO2', 'Ta2O5', 'Air'],
        'Thickness_List': [0, 1e-7, 1e-7, 0 ],
        'Lambda_List': [startwl, endwl, numwl],
        'Temperature': temp
       }
### create multilayer
dummy = multilayer(structure)
lambdas = dummy.lambda_array

### find Lamba_List indices for wl1 and wl2
d_lamb = (endwl - startwl)/(numwl-1)
wl1i = round((wl1 - startwl)/(d_lamb))
wl2i = round((wl2 - startwl)/(d_lamb))

### get refractive index arrays
n_SiO2 = dummy.layer_ri(1)
n_Ta2O5 = dummy.layer_ri(2)

### calculate effective quarter wave lengths for given vacuum wavelengths
thickL1 = 0.25*wl1/n_SiO2[wl1i].real
thickL2 = 0.25*wl2/n_SiO2[wl2i].real
thickH1 = 0.25*wl1/n_Ta2O5[wl1i].real
thickH2 = 0.25*wl2/n_Ta2O5[wl2i].real

print("thickL1: ",thickL1)
print("thickL2: ",thickL2)
print("thickH1: ",thickH1)
print("thickH2: ",thickH2)

### create a full filter-PS-filter stack
full_structure = {
        'Material_List' : ['Air', 'SiO2', 'Al', 'SiO2', 'Ta2O5', 'SiO2', 'Ta2O5', 'SiO2', 'PS', 'SiO2', 'Al', 'SiO2', 'Ta2O5', 'SiO2', 'Ta2O5', 'SiO2', 'Air'],
        'Thickness_List': [0, thickL1/2, 1.5e-9, thickL1/2, thickH1, thickL1, thickH1, 1e-4, 1e-3, thickL2/2, 1.5e-9, thickL2/2, thickH2, thickL2, thickH2, 1e-4, 0],
        'Lambda_List': [startwl, endwl, numwl],
        'Temperature': temp
       }

### create the bottom filter stack
structure_bottom = {
        'Material_List' : ['Air', 'SiO2', 'Al', 'SiO2', 'Ta2O5', 'SiO2', 'Ta2O5', 'SiO2', 'Air'],
        'Thickness_List': [0, thickL2/2, 1.5e-9, thickL2/2, thickH2, thickL2, thickH2, 1e-4, 0],
        'Lambda_List': [startwl, endwl, numwl],
        'Temperature': temp
       }
# create top filter
structure_top = {
        'Material_List' : ['Air', 'SiO2', 'Al', 'SiO2', 'Ta2O5', 'SiO2', 'Ta2O5', 'SiO2', 'Air'],
        'Thickness_List': [0, thickL1/2, 1.5e-9, thickL1/2, thickH1, thickL1, thickH1, 1e-4, 0],
        'Lambda_List': [startwl, endwl, numwl],
        'Temperature': temp
       }


### create multilayer
full = multilayer(full_structure)

### create multilayer
bottom = multilayer(structure_bottom)
top = multilayer(structure_top)



### Get the solar spectrum
AM = datalib.AM(full.lambda_array)

### Get the solar flux spectrum (number of photons per second per wavelength per meter squared)
AMflux = AM * full.lambda_array / (h * c)

### Integrate the above-gap solar flux transmitted through optic
AM_transmit = numlib.Integrate(AM * (1 - full.reflectivity_array), full.lambda_array, startwl, 760e-9)
print(AM_transmit," photons / s / m^2 transmitted")

### Integrate the total above-gap solar flux
AM_tot = numlib.Integrate(AM, full.lambda_array, startwl, 760e-9)
print(AM_tot," photons / s / m^2 total from Sun")

print(100* AM_transmit / AM_tot,"% transmissive")


### Make a BB spectrum at specified temp
rho = datalib.BB(full.lambda_array, temp)

### plot reflectivity and absorptivity/emissivity of the multilayer


### Integrate the blocked BB radiation longer than 3.5 um
BB_blocked = numlib.Integrate(rho * bottom.reflectivity_array, bottom.lambda_array, 3500e-9, endwl)
print(BB_blocked," photons / s / m^2 blocked")


### Integrate the total BB radiation longer than 3.5 um
BB_tot = numlib.Integrate(rho, full.lambda_array, 3500e-9, endwl)
print(BB_tot," photons / s / m^2 total from BB radiation")

print(100* BB_blocked / BB_tot,"% BB blocked")

### Integrate the blocked BB radiation longer than 3.5 um
BB_trapped = 0.5*numlib.Integrate(rho * bottom.reflectivity_array, bottom.lambda_array, 3500e-9, endwl) + 0.5*numlib.Integrate(rho * top.reflectivity_array, top.lambda_array, 3500e-9, endwl)
print(BB_trapped," photons / s / m^2 blocked")

print(100* BB_trapped / BB_tot,"% BB trapped")

### Integrate the transmitted BB radiation shorter than 3.5 um
EBE_transmitted = numlib.Integrate(rho * (1 - bottom.reflectivity_array), bottom.lambda_array, 1000e-9, 3500e-9)
print(EBE_transmitted," photons / s / m^2 above exciton biding energy transmitted")

### Integrate the total BB radiation shorter than 3.5 um
EBE_tot = numlib.Integrate(rho, full.lambda_array, 1000e-9, 3500e-9)
print(EBE_tot," photons / s / m^2 total above exciton binding energy from BB radiation")

print(100* EBE_transmitted / EBE_tot,"% above binding energy transmitted")
print(structure)

ref_tot = numlib.Integrate(full.reflectivity_array, full.lambda_array, 3500e-9, endwl)
print(" reflection total of full stack 1:", ref_tot)

ref_bottom = numlib.Integrate(bottom.reflectivity_array, bottom.lambda_array, 3500e-9, endwl)
ref_top = numlib.Integrate(top.reflectivity_array, top.lambda_array, 3500e-9, endwl)
print(" reflection total of bottom stack 1:", ref_bottom)
print(" reflection total of top stack 1:", ref_top)


print("structure 1:",full.d)
print("bottom 1:", bottom.d)
### create a full filter-PS-filter stack
full_structure = {
        'Material_List' : ['Air', 'SiO2', 'Al', 'SiO2', 'Ta2O5', 'SiO2', 'Ta2O5', 'SiO2', 'PS', 'SiO2', 'Al', 'SiO2', 'Ta2O5', 'SiO2', 'Ta2O5', 'SiO2', 'Air'],
        'Thickness_List': [0, thickL1/2, 1.5e-9, thickL1/2, thickH1, thickL1, thickH1,  1e-4, 1e-3, thickL2/2, 1.5e-9, thickL2/2, thickH2, thickL2, thickH2, 1e-4, 0],
        'Lambda_List': [startwl, endwl, numwl],
        'Temperature': temp
       }
thickL1 = 0.25*wl1/n_SiO2[wl1i].real
thickL2 = 0.25*wl2/n_SiO2[wl2i].real
thickH1 = 0.25*wl1/n_Ta2O5[wl1i].real
thickH2 = 0.25*wl2/n_Ta2O5[wl2i].real
thickL1 *= 1e9
thickL2 *= 1e9
thickH1 *= 1e9
thickH2 *= 1e9


def update_bottom(x, rho):
    for i in range(0,len(x)):
        if i<7 and i!=1 and i!=6:
            top.d[i+1] = x[i] * 1e-9
        if i!=1 and i!=6 and i!=7 and i!=9 and i!=14: 
            full.d[i+1] = x[i] * 1e-9
        #if i>7 and i!=9 and i!=14:
        #    bottom.d[i-6] = x[i] * 1e-9
    bottom.d[1:8] = full.d[9:16]
    ### now we have the new structure, update fresnel quantities
    bottom.fresnel()
    top.fresnel()
    full.fresnel()
    print("structure 2:",full.d)
    print("bottom 2:",bottom.d)

    ### now we have new emissivity, update thermal emission
    BB_blocked = numlib.Integrate(rho * bottom.reflectivity_array, bottom.lambda_array, 3500e-9, 20300e-9)
    ### Integrate the total BB radiation longer than 3.5 um
    BB_tot = numlib.Integrate(rho, full.lambda_array, 3500e-9, endwl)
    ref_tot = numlib.Integrate(full.reflectivity_array, full.lambda_array, 3500e-9, endwl)
    print(" reflection total of full stack 2:", ref_tot)
    ref_bottom = numlib.Integrate(bottom.reflectivity_array, bottom.lambda_array, 3500e-9, endwl)
    ref_top = numlib.Integrate(top.reflectivity_array, top.lambda_array, 3500e-9, endwl)
    print(" reflection total of top 2:", ref_top)
    print(" reflection total of bottom 2:", ref_bottom)

    BB_trapped = 0.5*numlib.Integrate(rho * bottom.reflectivity_array, bottom.lambda_array, 3500e-9, endwl) + 0.5*numlib.Integrate(rho * top.reflectivity_array, top.lambda_array, 3500e-9, endwl)
    ### Integrate the transmitted BB radiation shorter than 3.5 um
    EBE_transmitted = numlib.Integrate(rho * (1 - bottom.reflectivity_array), bottom.lambda_array, 1000e-9, 3500e-9)
    ### Integrate the total BB radiation shorter than 3.5 u
    EBE_tot = numlib.Integrate(rho, full.lambda_array, 1000e-9, 3500e-9)
    print("BB_blocked: ",BB_blocked/BB_tot * 100)
    print("BB_trapped: ",BB_trapped/BB_tot * 100)
    print("EBE_transmitted: ",EBE_transmitted/EBE_tot * 100)
    return -0.1 * BB_blocked/BB_tot - 0.1 * BB_trapped/BB_tot - 0.8 * EBE_transmitted/EBE_tot

def BuildGradient(x0, rho):
    dim = len(x0)
    h0 = 0.5*np.ones(dim)
    g = np.zeros(dim)
    for i in range(0,dim):
        if i!=1 and i!=6 and i!=7 and i!=9 and i!=14:
            xpass = np.copy(x0)
            fx = x0[i] + h0[i]
            bx = x0[i] - h0[i]
            xpass[i] = fx
            efx = update_bottom(xpass, rho)
            xpass[i] = bx
            ebx = update_bottom(xpass, rho)
            run = 2*h0[i]
            g[i] = (efx-ebx)/run
    return g

### Function that gets the negative of the efficiency and the 
### negative of the gradient for use in the l-bfgs-b algorithm
### also prints out the time for timing purposes!
def SuperFunc(x0):
    en = update_bottom(x0, rho)
    c_time = time.time()
    print(en,",",c_time)
    gr = BuildGradient(x0, rho)
    return en, gr



x = np.array([thickL1/2, 1.5, thickL1/2, thickH1, thickL1, thickH1, 1e5,  1e6, thickL2/2, 1.5, thickL2/2, thickH2, thickL2, thickH2, 1e5])

print(x)

print(update_bottom(x, rho))
#g = BuildGradient(x, rho)
#print(g)


# the bounds for L-BFGS-B updates!
bfgs_xmin = np.ones(len(x))
bfgs_xmax = 800*np.ones(len(x))

# rewrite the bounds in the way required by L-BFGS-B
bfgs_bounds = [(low, high) for low, high in zip(bfgs_xmin, bfgs_xmax)]

bfgs_bounds[1] = (1.4,1.5)
bfgs_bounds[6] = (9.9e4, 1e5)
bfgs_bounds[7] = (9.9e5, 1e6)
bfgs_bounds[9] = (1.4, 1.5)
bfgs_bounds[14] = (9.9e4, 1e5)

print(bfgs_bounds)

### initialize the solution vector xs to be the thicknesses from 
### the structure dictionary
#xs = np.zeros(length)
#for i in range(0,length):
#    xs[i] = cc.d[i+1]*1e9
#xs = np.copy(x)
xs = np.copy(x)
for i in range(0,len(xs)):
    if i!=1 and i!=6 and i!=7 and i!=9 and i!=14:
        xs[i] = np.random.randint(low=1, high=800, size=1)
        
print("random layer!",xs)

### print out initial solution vector and initial efficiency
#print("xs is ")
#print(xs)
print("efficiency is ",update_bottom(xs, rho))

### run l-bfgs-b algorithm!
ret = minimize(SuperFunc, xs, method="L-BFGS-B", jac=True, bounds=bfgs_bounds)

### print optimal solution and its efficiency!
print(ret.x)
#print(update_bottom(ret.x, rho))
