#!/usr/bin/env python
# coding: utf-8

# In[1]:


from wptherml.wpml import multilayer
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import basinhopping
import time
import numpy as np


structure = {

        'Material_List' : ['Air','J-Agg', 'TiO2', 'AlN','TiO2', 'AlN','Ag', 'Air'],
        ### Thicknesses just chosen arbitrarily, replace with "optimal" values
        'Thickness_List': [0, 15e-9, 8e-9, 8e-9, 8e-9, 8e-9, 300e-9, 0 ],
        'Lambda_List': [300e-9, 1500e-9, 1000],
        'Temperature': 300,
        'Gradient_List':[2,3,4,5],
        'LIGHTBULB': 1,

        }


### create instance of multilayer structure
cc = multilayer(structure)

### how many degrees of freedom will we vary over?
length = len(cc.gradient_list)


# Define a function *update_multilayer(x0)* that will take an array of thicknesses *x0* 
# (in nm) for the layers-to-be-varied, update the multilayer with those new
# thicknesses, and return the enhancement factor. Also define a function
# *analytic_grad(x0)* that works the same way, but returns the gradient with 
# respect to the layer thicknesses of the layers-to-be-varied.

# In[2]:


def update_multilayer(x0):
    dim = len(x0)
    ### update each layer-to-be-varied
    for i in range(2,dim+2):
        cc.d[i]= x0[i-2]*1e-9
    ### recompute fresnel quantities
    cc.fresnel()
    ### recompute enhancement factor
    cc.jagg_enhancement()
    ### return negative of enhancement factor 
    ### recall that scipy's optimize functions will find minimum, 
    ### so we need to give them the negative of the objective
    ### we wish to maximize
    return -cc.jagg_enhancement_val

def analytic_grad(x0):
    dim = len(x0)
    g = np.zeros(dim)
    ### update multilayer and fresnel quantities
    cur = update_multilayer(x0)
    ### update gradient of fresnel quantities
    cc.fresnel_prime()
    ### compute gradient of objective
    cc.jagg_enhancement_prime()
    g = cc.jagg_enhancement_grad
    
    return -g*1e-9

def numeric_grad(x0):
    dim = len(x0)
    h0 = 0.01*np.ones(dim)
    g = np.zeros(dim)
    for i in range(0,dim):
        xpass = np.copy(x0)
        fx = x0[i] + h0[i]
        bx = x0[i] - h0[i]
        xpass[i] = fx
        efx = update_multilayer(xpass)
        xpass[i] = bx
        ebx = update_multilayer(xpass)
        run = 2*h0[i]
        g[i] = (efx-ebx)/run
    return g

### function that calls both update_multilayer and analytic_gradient
### and returns both objective and gradient
def SuperFunc1(x0):
    en = update_multilayer(x0)
    gr = analytic_grad(x0)
    return en, gr

### function that calls both update_multilayer and numerical_gradient
### and returns both objective and gradient
def SuperFunc(x0):
    en = update_multilayer(x0)
    gr = numeric_grad(x0)
    return en, gr

### Just to confirm that both 
### analytic and numerical gradients give the
### same results to within acceptable precision
x0 = np.zeros(2)
### this is geometry that maximizes enhancement, 
### so gradient should nearly vanish
x0[0] = 7.152392
x0[1] = 21.5937
en1, gr1 = SuperFunc1(x0)
en2, gr2 = SuperFunc(x0)
print("Gr1",en1,gr1)
print("Gr2",en2,gr2)

### this is a random point, gradient probably won't vanish
x0[0] = 13.1
x0[1] = 42.1
en1, gr1 = SuperFunc1(x0)
en2, gr2 = SuperFunc(x0)
print("Gr1",en1,gr1)
print("Gr2",en2,gr2)


# We will customize a few things related to the Basin Hopping algorithm

# In[3]:


### This function will be called by the BH algorithm 
### each time a new local optimum is found... will
### print geometry, objective, and absolute time
def print_fun(x, f, accepted):
    c_time = time.time()
    ### note: print x[0]... x[L]
    print("!!!!!!!!!!!!!!  at minimum:",x, f, int(accepted),c_time)

### This defines how new initial configurations are chosen by the BH
### for each sub-optimization
def my_take_step(x):
    xnew = np.copy(x)
    dim = len(xnew)
    ### pick random numbers for position
    for i in range(0,dim):
        rn = 30*np.abs(np.random.randn())
        xnew[i] = rn
    return xnew

### specifically defines bounds for Basin Hopping
### algorithm... if a sub-optimization lands outside of
### these bounds and converges, the solution will not be accepted
class MyBounds(object):
      ### note xmax and xmin need to have as many elements as there are thicknesses that are varied
    def __init__(self, xmax=35*np.ones(length), xmin=np.ones(length)):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin

### the bounds for L-BFGS-B updates!  If an update 
### within a L-BFGS updates takes you outside these bounds,
### the update step will be modified to prevent going out of bounds
bfgs_xmin = np.ones(length)
bfgs_xmax = 400*np.ones(length)

# rewrite the bounds in the way required by L-BFGS-B
bfgs_bounds = [(low, high) for low, high in zip(bfgs_xmin, bfgs_xmax)]

minimizer_kwargs = {"method": "L-BFGS-B", "jac": True, "bounds": bfgs_bounds}
mybounds = MyBounds()

### note: initialize xs to be length L where L are the number of layers to be varied
xs = np.ones(length)*20
print("xs is ",xs)

ret = basinhopping(SuperFunc, xs, minimizer_kwargs=minimizer_kwargs, niter=60, take_step=my_take_step, callback=print_fun, accept_test=mybounds)

print(ret.x)
print(update_multilayer(ret.x))


# In[ ]:




