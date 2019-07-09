.. raw:: html

   <img src="../../Logo/WPtherml.png" alt="drawing" width="200"/> 
   Pioneering the design of materials for harnessing heat.

Overview
--------

WPTherml stands for **W**\ illiam **P**\ aterson University's tool for
**Th**\ ermal **E**\ nergy and **R**\ adiation management with
**M**\ ulti **L**\ ayer nanostructures. The vision of this software
package is to provide an easy-to-use platform for the design of
materials with tailored optical and thermal properties for the vast
number of energy applications where control of absorption and emission
of radiation, or conversion of heat to radiation or vice versa, is
paramount. The optical properties are treated within classical
electrodynamics, and the current version uses the Transfer Matrix Method
to rigorously solve Maxwell's equations for layered isotropic media.
WPTherml was conceived and developed by the `Foley Lab`_ at William
Paterson University. More details of the Transfer Matrix equations,
along will the full mathematical formulation currently implemented in
WPTherml, can be found in the `documentation`_.

Quick Start
-----------

-  WPTherml is written in Python 3 and requires the numpy, scipy, and
   matplotlib packages. Current installation of the Anaconda Python 3
   package should provide all you need on Windows, Mac, or Linux
   platforms
-  To get started, clone or download this repository to your computer
-  Open a new .py file in your favorite text editor or IDE, e.g.

``vim example.py``

The capabilities of this package are contained within a class called
multilayer. A basic example of a script that imports the multilayer
class, computes the reflectivity of 20 nm gold film coated with 50 nm of
TiO2 and 100 nm SiO2, and plots it using pyplot follows:

.. code:: python

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

.. _Foley Lab: https://foleylab.github.io
.. _documentation: https://github.com/FoleyLab/wptherml/blob/master/docs/Equations.pdf
