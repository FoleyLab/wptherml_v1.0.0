
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

q = 1.60217662e-19
c=299792458
h=6.626e-34
k=1.38064852e-23


### functions for refractive index of different common dielectrics
def Material_RI(lam, arg):
    ## array of wavelengths should be in meters
    ## but we might need it in nm or microns
    l_nm = lam*1e9
    lmic = lam*1e6
    if (arg=='HfO2'):
        A = 187178
        B = 9993.46
        C = 0.0173801
        D = 1.939999
        n = A/(l_nm**4) + B/(l_nm**2) + C/l_nm + D + 0j/l_nm
    elif (arg=='Al2O3'):
        A = 187178
        B = 9993.46
        C = 0.0173801
        D = 1.69999
        n = A/(l_nm**4) + B/(l_nm**2) + C/l_nm + D + 0j/l_nm
    ### This model works well for glass into near IR
    elif (arg=='SiO2' and lam[len(lam)-1]<5000e-9):
        A = 187178
        B = 9993.46
        C = 0.0173801
        D = 1.45
        n = A/(l_nm**4) + B/(l_nm**2) + C/l_nm + D + 0j/l_nm
    #elif (arg=='TiO2') and lam[len(lam)-1]<5000e-9:
    #    A = 187178
    #    B = 9993.46
    #    C = 0.0173801
    #    D = 2.4
    #    n = A/(l_nm**4) + B/(l_nm**2) + C/l_nm + D + 0j/l_nm
    elif (arg=='AlN' and lam[len(lam)-1]<10000e-9):
        A = 1.859
        B = 0.3401
        n = A + B/(lmic*lmic) + 0j/lmic
    elif (arg=='Air'):
        A = 0.
        n = A/lam + 0j/lam + 1
    elif (arg=='TiN'):
        n = TiN_Drude_Lorentz(lam)
    elif (arg=='W' or arg=='HfN' or arg=='Re' or arg=='Rh' or arg=='Ru'):
        n = Read_RI_from_File(lam, arg)
    elif (arg=='Ag' or arg=='Au' or arg=='Pd' or arg=='Pt' or arg=='SiO2'):
        n = Read_RI_from_File(lam, arg)
    elif (arg=='AlN' or arg=='Si' or arg=='TiO2'):
        n = Read_RI_from_File(lam, arg)
    elif (arg=='J-Agg'):
        n = TDBC(lam)
    ### default is air    
    else:
        A = 0.
        n = A/lam + 0j/lam + 1
    return n

def TiN_Drude_Lorentz(lam):
    ci = 0+1j
    epsinf = 3.59
    rho = 1.09e-6
    eps0 = 8.854e-12
    hbar = 6.582e-16
    tau = 7.82e-16
    amp = 4.658129
    br = 4.3827
    en = 5.778

    l_nm = lam*1e9
    E = 1240./l_nm
    eps = epsinf - hbar*hbar/(eps0*rho*(tau*E*E + ci*hbar*E))
    eps = eps + amp*br*en/(en*en - E*E - ci*E*br)
    return np.sqrt(eps)



def Read_RI_from_File(lam, matname):
    if (matname=='W'):
        a = np.loadtxt('wptherml/datalib/W_Palik_RI_f.txt')
    elif (matname=='TiO2'):  #Need to re-order HfN data
        a = np.loadtxt('wptherml/datalib/TiO2_Siefke.txt')
    elif (matname=='Re'):
        a = np.loadtxt('wptherml/datalib/Re_Palik_RI_f.txt')
    elif (matname=='Ru'):
        a = np.loadtxt('wptherml/datalib/Ru_Palik_RI_f.txt')
    elif (matname=='Rh'):
        a = np.loadtxt('wptherml/datalib/Rh_Palik_RI_f.txt')
    elif (matname=='Ag' and lam[len(lam)-1]<=1000e-9):
        a = np.loadtxt('wptherml/datalib/Ag_JC_RI_f.txt')
    elif (matname=='Ag' and lam[len(lam)-1]>1000e-9):
        a = np.loadtxt('wptherml/datalib/Ag_Yang.txt')
    elif (matname=='Au' and lam[len(lam)-1]<=1000e-9):
        a = np.loadtxt('wptherml/datalib/Au_JC_RI_f.txt')
    elif (matname=='Au' and lam[len(lam)-1]>1000e-9):
        a = np.loadtxt('wptherml/datalib/Au_IR.txt')
    elif (matname=='Pd'):
        a = np.loadtxt('wptherml/datalib/Pd_Palik_RI_f.txt')
    elif (matname=='Pt'):
        a = np.loadtxt('wptherml/datalib/Pt_Palik_RI_f.txt')
    elif (matname=='SiO2'):
        a = np.loadtxt('wptherml/datalib/SiO2_IR.txt')
    elif (matname=='AlN'):
        a = np.loadtxt('wptherml/datalib/AlN_IR.txt')
    elif (matname=='Si'):
        a = np.loadtxt('wptherml/datalib/Si_Schinke.txt')
    elif (matname=='W_Al2O3_Alloy'):
        a = np.loadtxt('wptherml/datalib/W_Al2O3_Alloy.txt')
    else:
        a = np.loadtxt('wptherml/datalib/W_Palik_RI_f.txt')
    ### now that we have read in the text, interpolate/extrapolate RI
    datlam = np.zeros(len(a))
    datn   = np.zeros(len(a))
    datk   = np.zeros(len(a))
    for i in range(0,len(a)):
        datlam[i]  = a[i][0]
        datn[i] = a[i][1]
        datk[i] = a[i][2]
        
    ### use linear interpolation/extrapolation
    order = 1
    ### form the interpolator/extrapolator object for datn
    sn = InterpolatedUnivariateSpline(datlam, datn, k=order)
    ### form the interpolator/extrapolator object for datk
    sk = InterpolatedUnivariateSpline(datlam, datk, k=order)
    ### compute the interpolated/extrapolated values for real part of RI
    yn = sn(lam)
    ### compute the interpolated/extrapolated values for imaginary part of R
    yk = sk(lam)
    ### for complex RI array for each value of lambda
    n = yn + 1j*yk
    return n

### returns interpolated/extrapolated EQE of monocrystaline silicon PV cells
### given an input array of wavelengths
def SR_Si(lam):
    ### values of lambda along which EQE is experimeintally known
    datlam = np.linspace(260e-9, 1310e-9, 22)
    ### experimental values of EQE for monocrystaline Si... 
    dateqe = 0.01*np.array([0., 0., 11.5, 23., 33., 37., 41., 45., 49., 52., 56., 60., 64., 62.5, 51., 35., 27.5, 20., 12.5, 7.5, 0., 0.])
    ### order of the spline
    order = 1
    ### form the interpolator/extrapolator object
    s = InterpolatedUnivariateSpline(datlam, dateqe, k=order)
    ### compute the interpolated/extrapolated values
    y = s(lam)
    
    return y
    


### returns interpolated/extrapolated EQE of InGaAsSb PV cells given an input
### array of wavelengths
def EQE_InGaAsSb(lam):
    ### values of lambda along which EQE is experimentally known
    datlam = np.linspace(1000.0e-9,3000.0e-9,21)
    ### values of EQE that were experimentally measured
    dateqe = 0.01*np.array([48., 50., 52.5, 54., 58., 59., 60., 62., 61., 62., 61.5, 59., 54., 22., 2., 0., 0., 0., 0., 0., 0.])
    
    ### use linear interpolation/extrapolation
    order = 1
    ### form the interpolator/extrapolator object
    s = InterpolatedUnivariateSpline(datlam, dateqe, k=order)
    ### compute the interpolated/extrapolated values
    y = s(lam)

    ### uncomment to plot experimental values
    '''
    plt.figure()
    plt.plot(datlam, dateqe, 'o')
    ### plot interpolated/extrapolated values
    plt.plot(lam, y, 'blue')
    '''
    ### return extrapolated EQE
    return y

def SR_InGaAsSb(lam):

    ### values of lambda along which EQE is experimentally known
    #datlam = np.linspace(1000.0e-9,3000.0e-9,21)
    datlam = 1e-9*np.array([200,320,1000,1100,1200,1300,1400,1500,
                            1600,1700,1800,1900,2000,2100,2200,2300,2400,
                            2500,2600,2700,2800,2900,3000])
    ### values of EQE that were experimentally measured
    dateqe = 0.01*np.array([0.,0.,48., 50., 52.5, 54., 58., 59., 60., 62., 61., 62., 61.5, 59., 54., 22., 2., 0., 0., 0., 0., 0., 0.])
    datsr = q*dateqe*datlam/(h*c)
    ### use linear interpolation/extrapolation
    order = 1
    ### form the interpolator/extrapolator object
    s = InterpolatedUnivariateSpline(datlam, datsr, k=order)
    ### compute the interpolated/extrapolated values
    y = s(lam)

    ### uncomment to plot experimental values
    '''
    plt.figure()
    plt.plot(datlam, datsr, 'o')
    ### plot interpolated/extrapolated values
    plt.plot(lam, y, 'blue')
    '''
    ### return extrapolated EQE
    return y

### spectral response for Silicon 
    

### returns interpolated/extrapolated spectral response (A/W) of GaSb PV cells
### given an input array of wavelenghts
def SR_GaSb(lam):
    ### values of lambda along which SR is experimentally known
    datlam = np.linspace(400.0e-9,2000.0e-9,17)
    ### values of SR that were experimentally measured
    datsr = 0.01*np.array([0., 0., 17.,34.,49.,58.,61. ,68., 73., 80., 85., 88., 87., 70.,  2., 0.,  0.])
    ### use linear interpolation/extrapolation
    order = 1
    ### form the interpolator/extrapolator object
    s = InterpolatedUnivariateSpline(datlam, datsr, k=order)
    ### compute the interpolated/extrapolated values
    y = s(lam)
    
    ### uncomment to plot experimental values
    '''
    plt.figure()
    plt.plot(datlam, datsr, 'o')
    plt.plot(lam, y, '-')
    '''
    return y

def BB(lam, T):
    ### speed of light in SI
    c = 299792458
    ### plancks constant in SI
    h = 6.62607004e-34
    ### boltzmanns constant in SI
    kb = 1.38064852e-23
    
    rho = np.zeros_like(lam)
    
    for i in range(0,len(rho)):
        rho[i] = 2*h*c*c/lam[i]**5 * 1./(np.exp(h*c/(lam[i]*kb*T))-1)
     
    #plt.figure()
    #plt.plot(lam, rho, '-')
    return rho



### Example for how to call the EQE_InGaAsSb function with a custom range of wavelengths!
#lam  = np.linspace(10e-9,3000e-9,1000)
#qe = SR_InGaAsSb(lam)
#SR = SR_GaSb(lam)

#B = BB(lam, 3000)



###AM 1.5
  
def AM(lam):  ###lam is x SI is y
    a = np.loadtxt('wptherml/datalib/scaled_AM_1_5.txt')
    x = np.zeros(len(a))
    y = np.zeros(len(a))
    ###  issue was just that the AM1.5 data had wavelength
    ###  in nanometers and we were assuming it was in meters!
    ###  now it is converted to meters at the time that it
    ###  is stored to the array called x
    for i in range(0,len(a)):
        x[i] = a[i][0]
        y[i] = a[i][1]
    datlam = x
    dateqe = y
    order = 1
    s = InterpolatedUnivariateSpline(datlam, dateqe, k=order)
    z = s(lam)
    #plt.plot(x,y,'blue')
    #plt.plot(lam, z, 'r--')
    #plt.show()
    return z

def ATData(lam):
    a = np.loadtxt('wptherml/datalib/ATrans.txt')
    x = np.zeros(len(a))
    y = np.zeros(len(a))

    for i in range(0,len(a)):
        x[i] = a[i][0]*1e-6
        y[i] = a[i][1]
    #plt.plot(x,y)
    datlam = x
    dateqe = y
    order = 1
    s = InterpolatedUnivariateSpline(datlam, dateqe, k=order)
    z = s(lam)
    #plt.plot(lam,z)
    #plt.show()
    return z
    

''' something wrong with reading this pl.txt text file!
def PhLum(lam):
    a = np.loadtxt('wptherml/datalib/pl.txt')
    x = np.zeros(len(a))
    y = np.zeros(len(a))
    ###  issue was just that the AM1.5 data had wavelength
    ###  in nanometers and we were assuming it was in meters!
    ###  now it is converted to meters at the time that it
    ###  is stored to the array called x
    for i in range(0,len(a)):
        x[i] = a[i][0]*1e-9
        y[i] = a[i][1]
    datlam = x
    datph = y
    order = 1
    s = InterpolatedUnivariateSpline(datlam, datph, k=order)
    z = s(lam)
    plt.plot(x,y,'blue')
    plt.plot(lam, z, 'r--')
    plt.show()
    return z
    
'''
### Given an array of lambda values, this function 
### will evaluate a Gaussian fit to photopic luminosity function
### and return an array of values from that fit.
def PhLum(lam):
    ## if changed change light lib functions

    ### gaussian parameters determined from fit to Photopic luminosity 
    ### function data, where raw data was accessed from here: 
    ### http://www.cvrl.org/database/data/lum/linCIE2008v2e_5.htm
    ### and fit was performed with gnuplot with the function f(x) = a*exp(-b*(x-c)**2)
    a               = 1.02433
    b               = 2.59462e+14
    c               = 5.60186e-07
    ### should be able to evaluate function at all values of array in one line
    z = a*np.exp(-b*(lam-c)**2)
    return z

### Read in CIE color matching functions from data file and interpolate/extrapolate
### on lam
def CIE(lam):
    a = np.loadtxt('wptherml/datalib/cie-cmf.txt')
    l = np.zeros(len(a))
    x = np.zeros(len(a))
    y = np.zeros(len(a))
    z = np.zeros(len(a))

    for i in range(0,len(a)):
        l[i] = a[i][0]*1e-9
        x[i] = a[i][1] 
        y[i] = a[i][2]
        z[i] = a[i][3]
        
    ## now interpolate over lam range using linear interpolation
    order = 1
    ## red response
    xint = InterpolatedUnivariateSpline(l, x, k=order)
    xbar = xint(lam)
    ## green response
    yint = InterpolatedUnivariateSpline(l, y, k=order)
    ybar = yint(lam)
    ## blue response
    zint = InterpolatedUnivariateSpline(l, z, k=order)
    zbar = zint(lam)

    ### store all the interpolated color matching functions 
    ### in a dictionary and return
    cie = {"xbar": xbar, 
         "ybar": ybar, 
         "zbar": zbar } 

    return cie

def TDBC(lam):
    ### values of lambda along which SR is experimentally known
    datlamk = 1e-9*np.array([380, 390, 400, 425, 450, 475, 500, 525, 550, 575, 582, 600, 625, 650, 675, 700, 710, 720])
    ### values of SR that were experimentally measured
    datk = np.array([0, 0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.01, 0.02, 0.08, 0.11, 0.09, 0.02, 0.01, 0.008, 0.004, 0, 0])
    ### use linear interpolation/extrapolation
    order = 1
    ### form the interpolator/extrapolator object
    sk = InterpolatedUnivariateSpline(datlamk, datk, k=order)
    ### compute the interpolated/extrapolated values
    yk = sk(lam)
    
    datlamn = 1e-9*np.array([380, 390, 400, 425, 450, 475, 500, 525, 550, 566, 575, 600, 625, 650, 675, 700, 710, 720])
    datn = 	np.array([1.550, 1.550, 1.545, 1.540, 1.535, 1.527, 1.520, 1.510, 1.490, 1.475, 1.480, 1.595, 1.580, 1.565, 1.555, 1.550, 1.550, 
                      1.550])
    
    sn = InterpolatedUnivariateSpline(datlamn, datn, k=order)
    yn = sn(lam)

    ### uncomment to plot experimental values
    
    #plt.figure()
    #plt.plot(datlamk, datk, 'o')
    #plt.plot(lam, yk, '-')
    #plt.plot(datlamn, datn, 'o')
    #plt.plot(lam, yn, '-')
    ci = 0+1j
    return yn + ci*yk

    
