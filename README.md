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
1.  Key:  Lambda_List - Data: a list of three floats that includes in order (i) shortest wavelength in meters, (ii) longest wavelength in meters, and (iii) total number of wavelengths where you would like the optical quantities to be evaluated
2.  Key:  Thickness_List - Data: a list of floats that specify the thickness in meters of each layer.  Note that the terminal layers (first and last) must have thickness of 0.
3.  Key:  Material_List - Data: a list of strings that specify the materials in each layer.  The following strings are currently recognized for the following supported materials:
   * Air - Air
   * SiO2 - Glass
   * HfO2 - Halfnium Oxide
   * Al2O3 - Aluminum Oxide
   * TiO2 - Titanium Oxide
   * AlN  - Aluminum Nitride
   * TiN - Titanium Nitride
   * Ag - Silver
   * Au - Gold
   * Pd - Palladium
   * Pt - Platinum
   * W - Tungsten
  
