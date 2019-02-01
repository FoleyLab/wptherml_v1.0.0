# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 13:02:37 2018

@author: varnerj
"""

from wptherml import TMM
from wptherml.datalib import datalib
from wptherml.tmmcore import tmm_core 
from wptherml import stpvlib
import numpy as np
import matplotlib.pyplot as plt
from wptherml.numlib import numlib


read_array = np.loadtxt("datalib/Trans9_5.txt")


for i in range(0,len(read_array)):
    lam[i] = read_array[i][0]

plt.plot(ATData(lam))
plt.show()






