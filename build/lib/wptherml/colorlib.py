from wptherml.numlib import numlib
from wptherml.datalib import datalib
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx

### Given incident spectrum, will
### classify a color as either being
### red, orange, yellow, green, blue, indigo, violet
### based on "distance" in actual RGB value from the 
### accepted "RGB" values of these colors
def classify_color(spec, lam):
    ### difference vector for each color
    ### define basis vectors for each color!
    colors = {
            ### these values coming ~ from https://www.rapidtables.com/web/color/RGB_Color.html
            'red': [255/255., 0, 0],
            'orange': [255/255., 128/255., 0],
            'yellow': [255/255., 255/255., 0],
            'green':  [0, 255/255., 0],
            'blue': [0, 0, 255/255.],
            'indigo': [75/255., 0, 130/255.],
            'violet': [255/255., 0, 255/255.]
    }
    ### this list can access the colors dictionary keys
    color_keys = ['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet']
    ### get RGB vector for current spectrum
    rgb = RGB_FromSpec(spec, lam)
    ### distance list
    distance_list = np.zeros(len(color_keys))
    for i in range(0, len(distance_list)):
        temp = colors[color_keys[i]]
        distance = np.sqrt((rgb[0]-temp[0])**2 + (rgb[1]-temp[1])**2 + (rgb[2]-temp[2])**2)
        distance_list[i] = distance
    color = color_keys[np.argmin(distance_list)]
    return color
    
    
    
    
    
    

def RGB_FromSpec(TE, lam):

    ### get color response functions 
    cie = datalib.CIE(lam)
    
    ### get X from TE spectrum 
    X = numlib.Integrate(TE*cie['xbar'], lam, 380e-9, 780e-9)

    ### get Y from TE spectrum
    Y = numlib.Integrate(TE*cie['ybar'], lam, 380e-9, 780e-9)
    
    ### get Z from TE spectrum
    Z = numlib.Integrate(TE*cie['zbar'], lam, 380e-9, 780e-9)
    
    ### get total magnitude
    tot = X+Y+Z
    
    ### get normalized values
    x = X/tot
    y = Y/tot
    z = Z/tot 
    ## should also be equal to z = 1 - x - y
    ### array of xr, xg, xb, xw, ..., zr, zg, zb, zw
    ### use hdtv standard
    xrgbw = [0.670, 0.210, 0.150, 0.3127]
    yrgbw = [0.330, 0.710, 0.060, 0.3291]
    zrgbw = []
    for i in range(0,len(xrgbw)):
        zrgbw.append(1. - xrgbw[i] - yrgbw[i])
    #print("zrgbw is ",zrgbw)
    
    ## rx = yg*zb - yb*zg
    rx = yrgbw[1]*zrgbw[2] - yrgbw[2]*zrgbw[1]
    ## ry = xb*zg - xg*zb
    ry = xrgbw[2]*zrgbw[1] - xrgbw[1]*zrgbw[2]    
    ## rz = (xg * yb) - (xb * yg)
    rz = xrgbw[1]*yrgbw[2] - xrgbw[2]*yrgbw[1]
    ## gx = (yb * zr) - (yr * zb)
    gx = yrgbw[2]*zrgbw[0] - yrgbw[0]*zrgbw[2]
    ## gy = (xr * zb) - (xb * zr)
    gy = xrgbw[0]*zrgbw[2] - xrgbw[2]*zrgbw[0]
    ## gz = (xb * yr) - (xr * yb)
    gz = xrgbw[2]*yrgbw[0] - xrgbw[0]*yrgbw[2]
    ## bx = (yr * zg) - (yg * zr)
    bx = yrgbw[0]*zrgbw[1] - yrgbw[1]*zrgbw[0]
    ## by = (xg * zr) - (xr * zg)
    by = xrgbw[1]*zrgbw[0] - xrgbw[0]*zrgbw[1]
    ## bz = (xr * yg) - (xg * yr)
    bz = xrgbw[0]*yrgbw[1] - xrgbw[1]*yrgbw[0]
    
    
    ## rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
    rw = (rx * xrgbw[3] + ry * yrgbw[3] + rz * zrgbw[3])/yrgbw[3] 
    ## gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
    gw = (gx * xrgbw[3] + gy * yrgbw[3] + gz * zrgbw[3])/yrgbw[3]
    ## bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;
    bw = (bx * xrgbw[3] + by * yrgbw[3] + bz * zrgbw[3])/yrgbw[3]

    ## /* xyz -> rgb matrix, correctly scaled to white. */    
    rx = rx / rw  
    ry = ry / rw  
    rz = rz / rw
    gx = gx / gw  
    gy = gy / gw  
    gz = gz / gw
    bx = bx / bw  
    by = by / bw  
    bz = bz / bw
    
    
    ## /* rgb of the desired point */
    r = (rx * x) + (ry * y) + (rz * z)
    g = (gx * x) + (gy * y) + (gz * z)
    b = (bx * x) + (by * y) + (bz * z)

    rgblist = []
    rgblist.append(r)
    rgblist.append(g)
    rgblist.append(b)

    # are there negative values?
    w = np.amin(rgblist)
    if w<0:
        rgblist[0] = rgblist[0] - w
        rgblist[1] = rgblist[1] - w
        rgblist[2] = rgblist[2] - w

    # scale things so that max has value of 1
    mag = np.amax(rgblist)
    

    rgblist[0] = rgblist[0]/mag
    rgblist[1] = rgblist[1]/mag
    rgblist[2] = rgblist[2]/mag
    
    #rgb = {'r': rgblist[0]/mag, 'g': rgblist[1]/mag, 'b': rgblist[2]/mag }
    
    return rgblist

def FalseColor_FromSpec(TE, lam):

    ### get color response functions 
    #cie = datalib.CIE((lam+1500e-9)*400e-9/2400e-9)
    SR = datalib.SR_InGaAsSb(lam)
    #plt.plot(1e9*lam, cie['xbar'], 'red', 1e9*lam, cie['ybar'], 'green', 1e9*lam, cie['zbar'], 'blue')
    #plt.show()

    ### This will be the shifted color cone response model for InGaAsSb
    ### get X from TE spectrum 
    #X = numlib.Integrate(TE*SR*cie['xbar'], lam, 100e-9, 4000e-9)

    ### get Y from TE spectrum
    #Y = numlib.Integrate(TE*SR*cie['ybar'], lam, 100e-9, 4000e-9)
    
    ### get Z from TE spectrum
    #Z = numlib.Integrate(TE*SR*cie['zbar'], lam, 100e-9, 4000e-9)

    ### this will be the step-function model for the response of 
    ### InGaAsSb PV cell
    ### get X from TE spectrum only - red is a penalty coming from sub-bg emission 
    X = numlib.Integrate(TE, lam, 2250e-9, 10000e-9)

    ### get Y from TE spectrum weighted by SR function... green is good!
    Y = numlib.Integrate(TE*SR, lam, 1900e-9, 2250e-9)
    
    ### get Z from TE spectrum weighted by SR function... blue is less good!
    Z = numlib.Integrate(TE*SR, lam, 400e-9, 1900e-9)
    
    ### get total magnitude
    tot = X+Y+Z
    
    ### get normalized values
    x = X/tot
    y = Y/tot
    z = Z/tot 
    ## should also be equal to z = 1 - x - y
    ### array of xr, xg, xb, xw, ..., zr, zg, zb, zw
    ### use hdtv standard
    xrgbw = [0.670, 0.210, 0.150, 0.3127]
    yrgbw = [0.330, 0.710, 0.060, 0.3291]
    zrgbw = []
    for i in range(0,len(xrgbw)):
        zrgbw.append(1. - xrgbw[i] - yrgbw[i])
    #print("zrgbw is ",zrgbw)
    
    ## rx = yg*zb - yb*zg
    rx = yrgbw[1]*zrgbw[2] - yrgbw[2]*zrgbw[1]
    ## ry = xb*zg - xg*zb
    ry = xrgbw[2]*zrgbw[1] - xrgbw[1]*zrgbw[2]    
    ## rz = (xg * yb) - (xb * yg)
    rz = xrgbw[1]*yrgbw[2] - xrgbw[2]*yrgbw[1]
    ## gx = (yb * zr) - (yr * zb)
    gx = yrgbw[2]*zrgbw[0] - yrgbw[0]*zrgbw[2]
    ## gy = (xr * zb) - (xb * zr)
    gy = xrgbw[0]*zrgbw[2] - xrgbw[2]*zrgbw[0]
    ## gz = (xb * yr) - (xr * yb)
    gz = xrgbw[2]*yrgbw[0] - xrgbw[0]*yrgbw[2]
    ## bx = (yr * zg) - (yg * zr)
    bx = yrgbw[0]*zrgbw[1] - yrgbw[1]*zrgbw[0]
    ## by = (xg * zr) - (xr * zg)
    by = xrgbw[1]*zrgbw[0] - xrgbw[0]*zrgbw[1]
    ## bz = (xr * yg) - (xg * yr)
    bz = xrgbw[0]*yrgbw[1] - xrgbw[1]*yrgbw[0]
    
    
    ## rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
    rw = (rx * xrgbw[3] + ry * yrgbw[3] + rz * zrgbw[3])/yrgbw[3] 
    ## gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
    gw = (gx * xrgbw[3] + gy * yrgbw[3] + gz * zrgbw[3])/yrgbw[3]
    ## bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;
    bw = (bx * xrgbw[3] + by * yrgbw[3] + bz * zrgbw[3])/yrgbw[3]

    ## /* xyz -> rgb matrix, correctly scaled to white. */    
    rx = rx / rw  
    ry = ry / rw  
    rz = rz / rw
    gx = gx / gw  
    gy = gy / gw  
    gz = gz / gw
    bx = bx / bw  
    by = by / bw  
    bz = bz / bw
    
    
    ## /* rgb of the desired point */
    r = (rx * x) + (ry * y) + (rz * z)
    g = (gx * x) + (gy * y) + (gz * z)
    b = (bx * x) + (by * y) + (bz * z)

    rgblist = []
    rgblist.append(r)
    rgblist.append(g)
    rgblist.append(b)

    # are there negative values?
    w = np.amin(rgblist)
    if w<0:
        rgblist[0] = rgblist[0] - w
        rgblist[1] = rgblist[1] - w
        rgblist[2] = rgblist[2] - w

    # scale things so that max has value of 1
    mag = np.amax(rgblist)
    

    rgblist[0] = rgblist[0]/mag
    rgblist[1] = rgblist[1]/mag
    rgblist[2] = rgblist[2]/mag
    
    #rgb = {'r': rgblist[0]/mag, 'g': rgblist[1]/mag, 'b': rgblist[2]/mag }
    
    
    return rgblist

def RenderColor(TE, lam, string):
    
    fig, ax = plt.subplots()
    # The grid of visible wavelengths corresponding to the grid of colour-matching
    # functions used by the ColourSystem instance.

    # Calculate the black body spectrum and the HTML hex RGB colour string
    cierbg = RGB_FromSpec(TE, lam)
    #cierbg = [1.,0.427,0.713]
    # Place and label a circle with the colour of a black body at temperature T
    x, y = 0., 0.
    circle = Circle(xy=(x, y*1.2), radius=0.4, fc=cierbg)
    ax.add_patch(circle)
    ax.annotate(string, xy=(x, y*1.2-0.5), va='center',
                ha='center', color=cierbg)

    # Set the limits and background colour; remove the ticks
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor('k')
    #ax.set_axis_bgcolor('k')
    # Make sure our circles are circular!
    ax.set_aspect("equal")
    plt.show()

    return 1

def RenderColor_Deuteranopia(TE, lam, T):

    fig, ax = plt.subplots()

    ### this is the normal array of RGB values
    cierbg = RGB_FromSpec(TE, lam)
    cierbg = [1.,0.427,0.713]
    ### dueteranopia R 
    rb = 0.625*cierbg[0] + 0.375*cierbg[1] 
    ### deuteranopia G 
    gb = 0.70*cierbg[0] + 0.3*cierbg[1]
    ### deuteranopia B
    bb = 0.3*cierbg[1] + 0.7*cierbg[2]

    ### ueteranopia R gets shifted by 0.56667 R
    cbrbg = [rb, gb, bb]
    scal = np.amax(cbrbg)
    cbrbg = cbrbg / scal
    x, y = 0., 0.
    circle = Circle(xy=(x, y*1.2), radius=0.4, fc=cbrbg)
    ax.add_patch(circle)
    ax.annotate('{:4d} K'.format(T), xy=(x, y*1.2-0.5), va='center',
                ha='center', color=cbrbg)

    # Set the limits and background colour; remove the ticks
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor('k')
    #ax.set_axis_bgcolor('k')
    # Make sure our circles are circular!
    ax.set_aspect("equal")
    plt.show()

    return 1
    

def RenderFalseColor(TE, lam, T):
    
    fig, ax = plt.subplots()
    # The grid of visible wavelengths corresponding to the grid of colour-matching
    # functions used by the ColourSystem instance.

    # Calculate the black body spectrum and the HTML hex RGB colour string
    cierbg = FalseColor_FromSpec(TE, lam)

    
    # Place and label a circle with the colour of a black body at temperature T
    x, y = 0., 0.
    circle = Circle(xy=(x, y*1.2), radius=0.4, fc=cierbg)
    ax.add_patch(circle)
    ax.annotate('{:4d} K'.format(T), xy=(x, y*1.2-0.5), va='center',
                ha='center', color=cierbg)

    # Set the limits and background colour; remove the ticks
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor('k')
    #ax.set_axis_bgcolor('k')
    # Make sure our circles are circular!
    ax.set_aspect("equal")
    plt.show()

    return 1
