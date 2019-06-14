#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:18:27 2019

@author: jay
"""





# First set up the figure, the axis, and the plot element we want to animate
'''
fig = plt.figure()
ax = plt.axes(xlim=(0, 300), ylim=(0, 10))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,
'''

from wptherml.wpml import multilayer
from wptherml.datalib import datalib
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib import animation

### Re-define structure:
structure = {
        ### No material actually called "Lorentz", so just use SiO2 as a placeholder
        ### that we can change later
        'Material_List' : ['Air', 'AlN', 'W', 'Air'],
        ### Thicknesses just chosen arbitrarily, replace with "optimal" values
        'Thickness_List': [0, 100e-9, 900e-9, 0],
        'Lambda_List': [300e-9, 4000e-9, 1000],
        'Temperature': 2700,
        'LIGHTBULB': 1

        }

#cc = multilayer(structure)
#print(cc.luminous_efficiency_val*100)

a = np.loadtxt('BH_data.txt')
### now that we have read in the text, interpolate/extrapolate RI
datd= np.zeros(len(a))
datn= np.zeros(len(a))
end = len(a)
for i in range(0,len(a)):
    datd[i]  = a[i][0]
    datn[i] = a[i][1]
### a function to take an array of values x_1, ..., x_8 as defined above,
### re-define the multi-layer structure according to its values, and
### then update the luminous efficiency accordingly:
def update_multilayer(x):
    ### recall x_1 is omega_p, x_2 is omega_0, x_3 is gamma
    #cc.layer_lorentz(1, omega_p, omega_0, gamma)
    cc.d[1] = x[0]*1e-9
    #print("Thickness is ",cc.d[1])
    ### now we have the new structure, update fresnel quantities
    cc.fresnel()
    ### now we have new emissivity, update thermal emission
    cc.thermal_emission()

    ### now we have new thermal emission, update luminous efficiency
    cc.luminous_efficiency()

    ### return luminous efficiency
    #print(x[0], cc.luminous_efficiency_val*100)
    return cc.luminous_efficiency_val*100

def analytic_grad(x0):
    dim = len(x0)
    g = np.zeros(dim)
    cur = update_multilayer(x0)
    for i in range(1,dim+1):
        cc.fresnel_prime(i)
        cc.luminous_efficiency_prime()
        g[i-1] = cc.luminous_efficiency_prime_val
    #print("g is ",-g*1e-7)
    return -g*1e-7


### define lorentz parameters

def SuperFunc(x):
    en = update_multilayer(x)
    gr = analytic_grad(x)
    return en, gr





def print_fun(x, f, accepted):
    print("!!!!!!!!!!!!!!  at minimum %.4f,   %.4f accepted %d" % (x[0], f, int(accepted)))

def my_take_step(x):
    xnew = np.copy(x)
    dim = len(xnew)
    for i in range(0,dim):
        rn = 180*np.abs(np.random.randn())
        xnew[i] = rn
    return xnew


class MyBounds(object):
    def __init__(self, xmax=[300], xmin=[1] ):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin

from scipy.optimize import minimize
from scipy.optimize import basinhopping
minimizer_kwargs = {"method": "BFGS", "jac": True}
mybounds = MyBounds()
xs = np.array([150])

trial = np.linspace(10.,300.,60)
n_vs_d = np.zeros(len(trial))
g_vs_d = np.zeros(len(trial))
idx = 0


for t in trial:
    xs[0] = t
    val = update_multilayer(xs)
    print("val is ",val)
    n_vs_d[idx] = val
#   g = analytic_grad(xs)
#   g_vs_d[idx] = g[0]
    idx = idx + 1

#plt.plot(trial, n_vs_d)
#plt.ylim(0,10)
#plt.show()

def _update_plot(i, fig, scat):
    #x = datd[i]
    #y = -1*datn[i]
    if (i%12==0):
        idx = int(i/12)
        
    scat.set_offsets(([datd[idx], datn[idx]])) #, [50, i], [100, i]))
    print('Frames: %d' %i)
    return scat,

fig = plt.figure()
plt.plot(trial, n_vs_d)
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


'''def animate(i):
    
    x = datd[i]
    y = -1*datn[i]
    #y = np.sin(2*np.pi*x/2)*np.cos(i/4.)
    line.set_data(x, y)
    return line,



anim = animation.FuncAnimation(fig, animate, init_func=init,
	                               frames=end, interval=1, blit=True)
#anim.save('PIB_EE3.mp4', fps=15, extra_args=['-vcodec', 'libx264'])
        
plt.show()
'''