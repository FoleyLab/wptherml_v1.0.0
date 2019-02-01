
### Numerically integrates real function stored in array y 
### evaluated at points stored
### in x between the ranges x[i] = a and x[f] = b
### just uses the rectangle rule
def Integrate(y, x, a, b):
    som = 0.
    keep = 1
    idx = 0
    ### if x array is going from large to small, reverse it!
    if (x[len(x)-1]-x[0]<0):
        xnew = x[::-1]
        x = xnew
        ynew = y[::-1]
        y = ynew
    for i in range(0,len(x)):
        ### are we outside of range?
        if x[i]>b:
            break
        ### dx for all elements except last
        elif (i<(len(x)-1)):
            dx = x[i+1]-x[i]
            fx = y[i]
            if x[i]>=a:
                som = som + fx*dx
        else:
            dx = x[i]-x[i-1]
            fx = y[i]
            if x[i]>=a:
                som = som + fx*dx

    return som
        
        
