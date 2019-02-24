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
1.  Key:  Lambda_List - Data: a list of three floats that includes in order (i) shortest wavelength in meters, (ii) longest wavelength in meters, and (iii) total number of wavelengths where you would like the optical quantities to be evaluated.  (Default is [400e-9,6000e-9,1000])
2.  Key:  Thickness_List - Data: a list of floats that specify the thickness in meters of each layer.  Note that the terminal layers (first and last) must have thickness of 0. (Default is [0, 900e-9, 0].)
3.  Key:  Material_List - Data: a list of strings that specify the materials in each layer (Default is ['Air', 'W', 'Air'].  
The following strings are currently recognized for the following supported materials:
   * Air - keyword for Air
   * SiO2 - keyword for Glass
   * HfO2 - keyword for Halfnium Oxide
   * Al2O3 - keyword for Aluminum Oxide
   * TiO2 - keyword for Titanium Oxide
   * AlN  - keyword for Aluminum Nitride
   * TiN - keyword for Titanium Nitride
   * Ag - keyword for Silver
   * Au - keyword for Gold
   * Pd - keyword for Palladium
   * Pt - keyword for Platinum
   * W - keyword for Tungsten
4.  Key: Temperature - Data: a float specifying the temperature of the multi-layer structure in Kelvin.  (Default is 300 K)
5.  Key: PV_Temperature - Data: a float specifying the temperature of a PV cell in a (S)TPV device in Kelvin.  (Default is 300 K).
6.  Key: Ambient_Temperature - Data: a float specifying the ambient temperature in Kelvin for radiative cooling applications. (Default is 300 K).
7.  Key: STPV_EMIT - Data: an int where '1' means compute properties associated with (S)TPV emitters. (Default is 0, do not compute these quantities).
8.  Key: STPV_ABS - Data: an int where '1' means compute properties associated with STPV/Concentrated Solar absorbers. (Default is 0).
9.  Key: COOLING - Data: an int where '1' means compute properties associated with radiative cooling. (Default is 0).
10.  Key: LIGHTBULB - Data: an int where '1' means compute properties associated with incandescent sources. (Default is 0).
11.  Key: COLOR - Data: an int where '1' means compute and display the ambient and thermal color of a structure. (Default is 0).
12.  Key: EXPLICIT_ANGLE - Data: an int where '1' means compute the optical properties and thermal emission at a range of angles and, when applicable, compute performance properties with explicit angular dependence.  (Default is 0, meaning most quantities will be computed assuming the emissivity does not depend upon angle.)
13.  Key: DEG - Data: an int that specifies the number of different angles that will be considered 
in the calculation of optical and thermal emission properties as a function of angle. (Default is 7, which has been observed to give reasonably good accuracy when all angular integrals are performed using Gauss-Legendre quadrature).

## Method and attribute list for multilayer class
Given the input parameters specified above, the multilayer class uses different methods to compute properties relevant for thermal applications, and those properties are store attributes
of the multilayer object.  The following is a list of methods of the multilayer class and their related attributes:

1.  
