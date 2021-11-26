#Dispersion Relation
#Solved by Newton-Raphson Method
#by Marco Quiroz
import numpy as np

#function disper is implemented
#Inputs:
#T: Period in seconds
#h: Depth in meters
#Output: Wavenumber (1/meters)
def disper(T,h,tol=0.000001):
    g = 9.81 #(m/s)
    xkh0 = ((2*np.pi/T)**2)*h/g
    coth = 1/np.tanh(xkh0**(3/4))
    xkh  = xkh0*(coth)**(2/3)
    while True:
        th = np.tanh(xkh)
        ch = np.cosh(xkh)
        f  = xkh0 - xkh*th
        fprima = -xkh/ch**2-th
        dxkh = -f/fprima
        if np.abs(dxkh/xkh) <= tol:
            xk = xkh/h
            break
        else:
            xkh = xkh + dxkh
    return xk


#Example
#input values
T = 10 #Period (sec)
h = 10 #Depth (m)

k = disper(10,10) #Wavenumber (1/m)
L = 2*np.pi/k     #Wavelength (m)

print('Wavelength = ',np.round(L,2),'(m)')
