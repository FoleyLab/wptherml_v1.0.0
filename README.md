<img src="Logo/WPtherml.png" alt="drawing" width="200"/> 
Pioneering the design of materials for harnessing heat.

## Overview
WPTherml stands for **W**illiam **P**aterson University's tool for **Th**ermal **E**nergy and **R**adiation management with **M**ulti **L**ayer nanostructures.
The vision of this software package is to provide an easy-to-use platform for the design of materials with tailored optical and thermal properties for
the vast number of energy applications where control of absorption and emission of radiation, or conversion of heat to radiation or vice versa, is paramount.
The optical properties are treated within classical electrodynamics, and the current version uses the Transfer Matrix Method to rigorously solve Maxwell's equations
for layered isotropic media.  WPTherml was conceived and developed by the [Foley Lab](https://foleylab.github.io) at William Paterson University. More details of the Transfer Matrix equations, along will the full mathematical formulation currently implemented in WPTherml, can be found in
the [documentation](https://github.com/FoleyLab/wptherml/blob/master/docs/Equations.pdf).

## Quick Start
- WPTherml is written in Python 3 and requires the numpy, scipy, and matplotlib packages.  Current installation of the Anaconda Python 3 package should provide all you need 
on Windows, Mac, or Linux platforms
- To get started, clone or download this repository to your computer
- Open a new .py file in your favorite text editor or IDE, e.g.

`vim test.py`

The capabilities of this package are contained within a class called multilayer.  A basic example 
of a script that imports the multilayer class, computes the reflectivity of 20 nm gold film coated with 50 nm of 
TiO2 and 100 nm SiO2, and plots
it using pyplot follows:
```python
from wptherml.wpml import multilayer
from matplotlib import pyplot as plt

### dictionary that stores basic properties 
### of the multilayer structure you want to simulate
structure = {
        ### actual materials the structure is made from... note terminal layers are air and
   	### top-side layer (layer upon which light is incident) is SiO2.
        ### Refractive index values are stored in the attribute self.n
        'Material_List': ['Air', 'SiO2', 'TiO2', 'Au', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 100e-9, 50e-9, 20e-9,  0],
         ### range of wavelengths optical properties will be calculated for
         ### values are stored in the array self.lam
        'Lambda_List': [400e-9, 800e-9, 1000]
        }

### create the instance called coated_au_film
coated_au_film = multilayer(structure)

### create a plot of the reflectivity of the coated au film - use red lines
### the wavelengths are stored in SI units so we will multiply by 1e9 to 
### plot them in nanometers
plt.plot(1e9*coated_au_film.lambda_array, coated_au_film.reflectivity_array, 'red')
plt.show()
```

- Save this script and run it either in the terminal as

`python test.py`

where test.py is the name of the file you created, or if you were doing this in an IDE, execute it within your IDE!

The schematic that illustrates the above example is shown in the figure below. Note the ordering of the 
layers in the picture and how they are specified through Material_List and Thickness_List relative to 
the incident, reflected, transmitted, and thermally-emitted light.

<img src="docs/Convention.png" alt="drawing" width="500"/>


There are illustrative examples of using the features of the multilayer class contained in Jupyter notebooks within this repository, including:

- [Validation of Basic Optical Properties](https://github.com/FoleyLab/wptherml/blob/master/Validate_Fresnel.ipynb)

- [Examples of Computing Basic Optical Properties](https://github.com/FoleyLab/wptherml/blob/master/Example1.ipynb)

- [Modeling Incandescent Sources](https://github.com/FoleyLab/wptherml/blob/master/Example2.ipynb)

- [Modeling Radiative Cooling Surfaces](https://github.com/FoleyLab/wptherml/blob/master/Validate_Cooling.ipynb)

- [Video Demo for Radiative Cooling](https://youtu.be/LC4TrnB8JK4)

more will be added in the near future!


## Playlist
The developers of WPTherml compiled a thematic [Spotify Playlist called "Everything Thermal"](https://open.spotify.com/playlist/1Vb7MV4WwjOMMHLbrX4TNN); we hope it will inspire you to imagine new possibilities for 
harnessing heat and thermal radiation!

## Features List
1. Computes Reflectivity, Transmissivity, and Absorptivity/Emissivity spectrum of arbitrary multi-layered planar structures using the Transfer Matrix Method
2. Computes Thermal Emission spectrum at a given temperature of multi-layer structure as emissivity * Blackbody spectrum 
3. Computes solar power absorbed from absorptivity * AM1.5 spectrum
4. From the quantities above, the following performance-related quantities can be computed for various thermal-related applications:
   * Spectral Efficiency of (S)TPV Emitters for a given PV
   * Useful Power Density (S)TPV Emitters for a given PV
   * Short Circuit Current Density (S)TPV Emitter for a given PV
   * TPV Efficiency (S)TPV Emitters for a given PV
   * Absorber Efficiency for STPV Absorbers for a given concentration factor
   * Luminous Efficiency/Luminous Efficacy of Incandescent bulb filaments
   * Cooling Power for day-time radiative cooling for a given ambient temperature and temperature of the multi-layer
5. From optical quantities, the following analysis can be performed
   * Identify Surface Plasmon Polariton modes
   * Identify Perfectly Absorbing modes
   * Rendering of color of a multi-layer at cool temperatures and at elevated temperatures

The calculations of the quantities above are facilitated by a class called *multilayer*.  The *multilayer* class parses a dictionary for key 
structural data like the material and thicknesses that comprise the multi-layer structure being modeled, the types of applications one wants to
consider the multi-layer structure for.  The following is the complete list of dictionary keys the *multilayer* class will recognize, along with
the data the user can supply in association with each key:
```python
'Lambda_List' # a list of three floats that includes in order (i) shortest wavelength in meters, (ii) longest wavelength in meters, and (iii) total number of wavelengths where you would like the optical quantities to be evaluated.  (Default is [400e-9,6000e-9,1000])

'Thickness_List' # a list of floats that specify the thickness in meters of each layer.  Note that the terminal layers (first and last) must have thickness of 0. (Default is [0, 900e-9, 0].)

'Material_List' # a list of strings that specify the materials in each layer (Default is ['Air', 'W', 'Air'].  
The following strings are currently recognized for the following supported materials:
   * 'Air' - keyword for Air
   * 'SiO2' - keyword for Glass
   * 'HfO2' - keyword for Hafnium Oxide
   * 'Al2O3' - keyword for Aluminum Oxide
   * 'TiO2' - keyword for Titanium Oxide
   * 'AlN'  - keyword for Aluminum Nitride
   * 'TiN' - keyword for Titanium Nitride
   * 'Ag' - keyword for Silver
   * 'Au' - keyword for Gold
   * 'Pd' - keyword for Palladium
   * 'Pt' - keyword for Platinum
   * 'W' - keyword for Tungsten

'Temperature'  # a float specifying the temperature of the multi-layer structure in Kelvin.  (Default is 300 K)

'PV_Temperature' # a float specifying the temperature of a PV cell in a (S)TPV device in Kelvin.  (Default is 300 K).

'Ambient_Temperature' # a float specifying the ambient temperature in Kelvin for radiative cooling applications. (Default is 300 K).

'STPV_EMIT' # an int where '1' means compute properties associated with (S)TPV emitters. (Default is 0, do not compute these quantities).

'STPV_ABS' # an int where '1' means compute properties associated with STPV/Concentrated Solar absorbers. (Default is 0).

'COOLING' # an int where '1' means compute properties associated with radiative cooling. (Default is 0).

'LIGHTBULB' # an int where '1' means compute properties associated with incandescent sources. (Default is 0).

'COLOR' # an int where '1' means compute and display the ambient and thermal color of a structure. (Default is 0).

'EXPLICIT_ANGLE' # an int where '1' means compute the optical properties and thermal emission at a range of angles and, when applicable, compute performance properties with explicit angular dependence.  (Default is 0, meaning most quantities will be computed assuming the emissivity does not depend upon angle.)

'DEG' # an int that specifies the number of different angles that will be considered 
in the calculation of optical and thermal emission properties as a function of angle. (Default is 7, which has been observed to give reasonably good accuracy when all angular integrals are performed using Gauss-Legendre quadrature).
```

## Method and attribute list for multilayer class
Given the input parameters specified above, the *multilayer* class uses different methods to compute properties relevant for thermal applications, and those properties are stored as attributes
of the *multilayer* object.  The following is a list of methods of the *multilayer* class and their related attributes:


```python
	def inline_structure(structure):
       	### a method to parse input parameters from a dictionary (here called structure, all currently-supported dictionary 
       	### keys are defined above.  This method is called by the __init__ and defines the following attributes:

	self.lambda_array 	# the list of wavelengths in meters that will be used to evaluate optical and thermal spectra
	self.d		  	# the list of thicknesses that define the geometry of the multilayer
	self.matlist      	# the list of strings that specify the materials
	self.n		  	# the 2D arrays of refractive index values for each material for each wavelength (inner index specifies material, outter index wavelength)
	self.T_ml         	# the temperature of the multi-layer in Kelvin
	self.T_cell       	# the temperature of the PV cell in Kelvin
	self.T_amb      	# the ambient temperature in Kelvin
	self.stpv_emitter_calc  # the flag that determines if (S)TPV emitter properties will be computed
	self.stpv_absorber_calc # the flag that determines if (S)TPV absorber properties will be computed
	self.cooling_calc    	# the flag that determines if radiative cooling properties will be computed
	self.lightbulb_calc     # the flag that determines if incandescent properties will be computed
	self.color_calc 	# the flag that determines if colors will be rendered
	self.explicit_angle 	# the flag that determines if explicit angle-dependence of optical properties will be considered
	self.deg		# the number of different angles that will be computed for angle-dependent optical properties
  ```
In addition to the attributes that are explicitly set by parsing user input, several more attributes that are arrays will be 
allocated based on attributes defined by inline_structure:
```python
	### The following are always created
	self.reflectivity_array 	# initialized as an array of zeros the same length as self.lambda_array
	self.transmissivity_array	# initialized as an array of zeros the same length as self.lambda_array
	self.emissivity_array		# initialized as an array of zeros the same length as self.lambda_array
	self.thermal_emission_array	# initialized as an array of zeros the same length as self.lambda_array

	### The following are created if self.explicit_angle == 1
	self.x				# points from Gauss-Legendre grid of degree self.deg from 0 to 1
	self.t				# self.deg angles on Gauss-Legendre grid transformed to be between 0 and pi/2
	self.w				# self.deg weights from Gauss-Legendre grid transformed to be between 0 and pi/2

	self.reflectivity_array_p       # initialized as a 2D array of zeros, inner dimension same as self.deg and outter same as self.lambda_array
        self.reflectivity_array_s       # initialized as a 2D array of zeros, inner dimension same as self.deg and outter same as self.lambda_array
        self.transmissivity_array_p     # initialized as a 2D array of zeros, inner dimension same as self.deg and outter same as self.lambda_array
        self.transmissivity_array_s     # initialized as a 2D array of zeros, inner dimension same as self.deg and outter same as self.lambda_array
        self.emissivity_array_p         # initialized as a 2D array of zeros, inner dimension same as self.deg and outter same as self.lambda_array
        self.emissivity_array_s         # initialized as a 2D array of zeros, inner dimension same as self.deg and outter same as self.lambda_array
        self.thermal_emission_array_p   # initialized as a 2D array of zeros, inner dimension same as self.deg and outter same as self.lambda_array
        self.thermal_emission_array_s   # initialized as a 2D array of zeros, inner dimension same as self.deg and outter same as self.lambda_array
```

```python
''' Method to compute optical properties of reflectivity, transmissivity, and 
emissivity of structure as a function of wavelength assuming normal incidence '''
def fresnel()

### Upon execution, the following arrays are filled with their respective values
### for every wavelength in self.lambda_array
self.reflectivity_array
self.transmissivity_array
self.emissivity_array
```

```python
''' Method to compute optical properties of reflectivity, transmissivity, and 
emissivity of structure as a function of wavelength and angle, both p- and s-polarizations
are considered '''
def fresnel_ea()

### Upon execution, the following arrays are filled with their respective values
### for every wavelength in self.lambda_array and every angle in self.t
self.reflectivity_array_p
self.reflectivity_array_s
self.transmissivity_array_p
self.transmissivity_array_s
self.emissivity_array_p
self.emissivity_array_s
```

```python
''' Method to compute thermal emission spectrum of a structure at a given temperature;
note temperature specified by self.T_ml '''
def thermal_emission()

### Upon execution, the following arrays are computed for every wavelength in self.lambda_array
### for temperature given by self.T_ml
self.BBs   # Blackbody spectrum
self.thermal_emission_array ## thermal emission of structure defined as Blackbody * emissivity
```


```python
''' Method to compute thermal emission spectrum of a structure at a given temperature for a range of angles '''
def thermal_emission_ea()

### Upon execution, the following arrays are computed for every wavelength in self.lambda_array
### and every angle in self.t for temperature given by self.T_ml
self.thermal_emission_array_p ## thermal emission of structure defined as Blackbody * p-polarized emissivity
self.thermal_emission_array_s ## thermal emission of structure defined as Blackbody * s-polarized emissivity
```
```python
''' Method to compute optical properties of reflectivity, transmissivity, 
and emissivity as a function of angle for a given polarization self.pol and wavelength lambda_0 '''
def angular_fresnel(self, lambda_0)

### Upon execution, the following arrays are computed for 180 angles between 0 and pi/2
self.r_vs_theta # reflectivity
self.t_vs_theta # transmissivity
self.eps_vs_theta # emissivity
```

```python
''' The following three methods compute figures of merit relevant for STPV emitters for a given
    temperature self.T_ml, PV type self.PV and bandgap self.lbg, and PV temperature self.T_cell.
    These methods assume the emissivity does not change with angle, and perform an analytic
    integration over solid angles that make the computations much quicker, though also less realistic.'''
self.stpv_se() # compute the spectral efficiency and stores it in the attribute self.spectral_efficiency_val
self.stpv_pd() # computes the useful power density and stores it in the attribute self.power_density_val
self.stpv_etatpv() # computes the TPV emitter efficiency and stores it in the attribute self.tpv_efficiency_val
```

```python
''' The following methods compute figures of merit relevant for STPV emitters for a given
    temperature self.T_ml, PV type self.PV and bandgap self.lbg, and PV temperature self.T_cell.
    These methods explicitly account for the angular dependence of the emissivity, making these calculations
    more realistic but also more time consuming. '''
self.stpv_se_ea() # compute the spectral efficiency and stores it in the attribute self.spectral_efficiency_val
self.stpv_pd_ea() # computes the useful power density and stores it in the attribute self.power_density_val
self.stpv_etatpv_ea() # computes the TPV emitter efficiency and stores it in the attribute self.tpv_efficiency_val
```
```python
''' The following methods compute the absorber efficiency of a STPV or concentrated solar absorber at a 
    given temperature self.T_ml '''
def stpv_etaabs_ea() # computes absorber efficiency and stores it in the attribute self.absorber_efficiency_val
```
```python
''' method to render color of a structure from its thermal emission at a given temperature self.T_ml '''
def thermal_color()
''' method to render color of a structure from its reflection spectrum '''
def ambient_color()
''' method to render color in a +/- 5nm band around the wavelength lambda '''
def pure_color(lambda)
```

```python
''' Method to compute the luminous efficiency of a structure at temperature self.T_ml.
    Stores value to self.luminous_efficiency_val '''
def luminous_efficiency()

''' Method to compute the radiative cooling power of a structure at temperature self.T_ml in ambient
    temperature self.T_amb while being illuminated by the AM1.5 spectrum.  Upon execution, the relevant
    values are stored to the attributes self.radiative_power_val (this is the flux that cools the structure),
    self.atmospheric_power_val (part of flux that warms the structure) and self.solar_power_val (part of the flux 
    that warms the structure).'''
def cooling_power()


''' Method  to add a layer to the structure; material of the layer to be added will be specified by 'material' argument
    and thickness of the layer will be specified by the 'thickness' argument.  The layer will be inserted after
    the 'layer_number' layer.  The method will also update spectral and performance quantities after the layer is
    added; the instance name will be preserved after execution, so this is like a mutation operation.'''
def insert_layer(layer_number, material, thickness)

''' Method to extract the array of refractive index values associated with a specific layer; the method returns 
    this array.  '''
def layer_ri(layer_number)

''' Method to define the refractive index of an existing layer (specified by layer_number) as an alloy
    of material_1 and material_2 with a specified volume_fraction of material_1 in material_2 according
    to either the Maxwell-Garnett or the Bruggeman effective medium theory.  Using 'Bruggeman' as the
    argument for model will use Bruggeman's effective medium theory, while any other string will default
    to Maxwell-Garnett theory. Optical properties and performance figures are NOT updated upon execution of this method.'''
def layer_alloy(layer_number, volume_fraction, material_1, material_2, model)

''' Method to define the refractive index of an existing layer (specified by layer number) as a single
    complex number (specified by refractive_index_value) for all wavelengths.  Optical properties and performance figures are NOT updated upon execution of this method.'''
def layer_static_ri(layer_number, refractive_index_value)

''' Method to compute complex wavevector magnitude associated with the surface plasmon polariton mode on a given multi-layer
    structure at a wavelength specified by the int wavelength_index, where self.lambda_array[wavelength_index] returns
    the wavelength you are interested in in meters.  Upon completion, the spp wavevector is stored in
    self.spp_resonance_val '''
def find_spp(wavelength_index)

''' Method to compute complex wavevector magnitude associated with the perfectly absorbing mode on a given multi-layer
    structure at a wavelength specified by the int wavelength_index, where self.lambda_array[wavelength_index] returns
    the wavelength you are interested in in meters.  Upon completion, the pa wavevector is stored in
    self.pa_resonance_val '''
def find_pa()
```

