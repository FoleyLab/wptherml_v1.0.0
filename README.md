<img src="Logo/WPtherml.png" alt="drawing" width="200"/> 
Pioneering the design of materials for harnessing heat.


## Features List
1. Computes Transfer Matrix given the following inputs from the user given the following input:
   * Number of layers in a multi-layer
   * Thickness of each layer
   * Refractive index of each layer (complex quantity)
   * Frequency/wavelength/wavenumber at which the TMM should be evaluated
   * Angle of incidence at which TMM should be evaluated (complex quantity)
2. From the Transfer Matrix, the following optical quantities can be computed:
   * reflection amplitude (complex quantity)
   * transmission amplitude (complex quantity)
   * Far-field Reflection 
   * Far-field Transmission
   * Absorbance/emittance
   * Stored Energy
   * Electric Field in each layer (complex quantity)
3. From these optical properties, the following quantities for thermal applications can be computed:
   * Thermal emission spectrum at a given T
   * Solar Absorption spectrum
   * Net solar absorptivity at a given T
4. From the quantities for thermal applications, the following performance-related quantities can be computed:
   * Spectral Efficiency ((S)TPV Emitters; numerical integration required, Response Function of PV required)
   * Useful Power Density ((S)TPV Emitters; numerical integration required, Response Function of PV required)
   * Short Circuit Current Density ((S)TPV Emitter; numerical integration required, Response Function of PV required)
   * TPV Efficiency ((S)TPV Emitters; numerical integration required, Response Function of PV required, PV Cell temperature required)
   * Absorber Efficiency (STPV Absorbers; numerical integration required, Concentration Factor, AM1.5 Spectrum required, and Temperature of Absorber Required)
   * Luminous Efficiency/Luminous Efficacy (Incandescent bulb filaments; numerical integration required, photopic luminosity function required, temperature of filament required)
   * Cooling Power (Passive cooling applications, see https://www.nature.com/articles/nature13883.pdf for more details)
4. From optical quantities, can perform the following analysis:
   * Identify Surface Plasmon Polariton modes
   * Identify Perfectly Absorbing modes
   
#TO DO:
- Upload TMM Module Code to this repository (James)
- Find a data set for the AM 1.5 spectrum and upload it to this repository (James)
- Find a data set for the Photopic luminosity function and upload it to this repository (Reem)
- Identify scipy or numpy cubic spline interpolation functions, provide simple example in this repository (Derek)
- Find data sets for response functions of the following PV cells: GaSb, InGaSb, InGaAsSb (JJF) 
- Find data set for transmissivity of atmosphere (Noor)
- Identify formula to calculate Rate of Cooling of a structure for radiative cooling application (Noor)
- Create a conceptual roadmap for the software package! (JJF)
   
  
