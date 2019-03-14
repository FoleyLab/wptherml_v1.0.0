from wptherml import tmm
from wptherml import colorlib
from wptherml import coolinglib
from wptherml import stpvlib
from wptherml import lightlib
from wptherml.datalib import datalib
from matplotlib import pyplot as plt
from scipy import integrate

import numpy as np

lam = np.linspace(0,2500*1e-9,1000)

InGaAsSb = datalib.SR_InGaAsSb(lam)

GaSb = datalib.SR_GaSb(lam)

plt.xlabel('Wavelength' r' $(nm)$')
plt.ylabel('Spectral Response' r' $(A/W)$')

plt.plot(lam*1e9, GaSb, 'orange',label = 'GaSb', linewidth = 2.25)
plt.plot(lam*1e9, InGaAsSb,label = 'InGaAsSb', linewidth = 2.25)
plt.legend(loc = 'best')
plt.show()
