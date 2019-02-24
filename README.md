<img src="Logo/WPtherml.png" alt="drawing" width="200"/> 
Pioneering the design of materials for harnessing heat.


## Features List
1. Computes Reflectivity, Transmissivity, and Absorptivity/Emissivity spectrum of arbitrary multi-layered planar structures using the Transfer Matrix Method
2. Computes Thermal Emission spectrum at a given temperature of multi-layer structure as emissivity * Blackbody spectrum 
3. Computes solar power absorbed from absorptivity * AM1.5 spectrum
4. From the quantities above, the following performance-related quantities can be computed for various thermal-related applications:
   * Spectral Efficiency of ((S)TPV Emitters for a given PV
   * Useful Power Density ((S)TPV Emitters for a given PV
   * Short Circuit Current Density ((S)TPV Emitter for a given PV
   * TPV Efficiency ((S)TPV Emitters for a given PV
   * Absorber Efficiency for STPV Absorbers for a given concentration factor
   * Luminous Efficiency/Luminous Efficacy of Incandescent bulb filaments
   * Cooling Power for day-time radiative cooling for a given ambient temperature and temperature of the multi-layer
5. From optical quantities, the following analysis can be performed
   * Identify Surface Plasmon Polariton modes
   * Identify Perfectly Absorbing modes
   * Rendering of color of a multi-layer at cool temperatures and at elevated temperatures

The calculations of the quantities above are facilitated by a class called $multilayer$.  The multilayer class parses a dictionary for key 
structural data like the material and thicknesses that comprise the multi-layer structure being modeled, the types of applications one wants to
consider the multi-layer structure for.  The following is the complete list of dictionary keys the multilayer class will recognize, along with
the data the user can supply in association with each key:
```python
'Lambda_List' # a list of three floats that includes in order (i) shortest wavelength in meters, (ii) longest wavelength in meters, and (iii) total number of wavelengths where you would like the optical quantities to be evaluated.  (Default is [400e-9,6000e-9,1000])

'Thickness_List' # a list of floats that specify the thickness in meters of each layer.  Note that the terminal layers (first and last) must have thickness of 0. (Default is [0, 900e-9, 0].)

'Material_List' # a list of strings that specify the materials in each layer (Default is ['Air', 'W', 'Air'].  
The following strings are currently recognized for the following supported materials:
   * 'Air' - keyword for Air
   * 'SiO2' - keyword for Glass
   * 'HfO2' - keyword for Halfnium Oxide
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
Given the input parameters specified above, the multilayer class uses different methods to compute properties relevant for thermal applications, and those properties are store attributes
of the multilayer object.  The following is a list of methods of the multilayer class and their related attributes:


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
        def inline_structure(structure):
        ### a method to parse input parameters from a dictionary (here called structure, all currently-supported dictionary
        ### keys are defined above.  This method is called by the __init__ and defines the following attributes:

        self.lambda_array       # the list of wavelengths in meters that will be used to evaluate optical and thermal spectra
        self.d                  # the list of thicknesses that define the geometry of the multilayer
        self.matlist            # the list of strings that specify the materials
        self.n                  # the 2D arrays of refractive index values for each material for each wavelength (inner index specifies material, outter index wavelength)
        self.T_ml               # the temperature of the multi-layer in Kelvin
        self.T_cell             # the temperature of the PV cell in Kelvin
        self.T_amb              # the ambient temperature in Kelvin
        self.stpv_emitter_calc  # the flag that determines if (S)TPV emitter properties will be computed
        self.stpv_absorber_calc # the flag that determines if (S)TPV absorber properties will be computed
        self.cooling_calc       # the flag that determines if radiative cooling properties will be computed
        self.lightbulb_calc     # the flag that determines if incandescent properties will be computed
        self.color_calc         # the flag that determines if colors will be rendered
        self.explicit_angle     # the flag that determines if explicit angle-dependence of optical properties will be considered
      
  
