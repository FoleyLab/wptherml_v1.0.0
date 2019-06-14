#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 10:34:58 2019

@author: jay
"""

from matplotlib import pyplot as plt
import numpy as np
from matplotlib import animation

a = np.loadtxt('BH_data.txt')
### now that we have read in the text, interpolate/extrapolate RI
datd= np.zeros(len(a))
datn= np.zeros(len(a))
end = len(a)
for i in range(0,len(a)):
    datd[i]  = a[i][0]
    datn[i] = a[i][1]

def _update_plot(i, fig, scat):
    #x = datd[i]
    #y = -1*datn[i]
    if (i%12==0):
        idx = int(i/12)
        
    scat.set_offsets(([datd[idx], datn[idx]])) #, [50, i], [100, i]))
    print('Frames: %d' %i)
    return scat,

fig = plt.figure()
#x = [0, 50, 100]
#y = [0, 0, 0]
x = [0]
y = [0]
ax = fig.add_subplot(111)
ax.grid(True,linestyle = '-', color = '0.75')
ax.set_xlim([0, 300])
ax.set_ylim([0,10])
scat = plt.scatter(x, y, c=x)
scat.set_alpha(0.8)

anim = animation.FuncAnimation(fig, _update_plot, fargs = (fig, scat), frames = 12*end, interval = 1)