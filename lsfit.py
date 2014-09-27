#!/usr/bin/env python
"""
google(Least squaries polynomial fitting in Python)
using numpy.linalg.lstsq
"""
import numpy as np
from numpy import pi,arctan,exp
import scipy as sp
import matplotlib.pyplot as plt

degree = 6

# generate a noisy parabola
t = np.linspace(0,100,200)
parabola = t**2
noise = np.random.normal(0,300,200)
y = parabola + noise

# form the Vandermonde matrix
A = np.vander(t,degree)

# find the x that minimize the norm of Ax-y
(coeffs, residuals, rank,sing_vals) = np.linalg.lstsq(A,y)

# create a polynomial using coefficients
f = np.poly1d(coeffs)

# for plot, estimate y for each observation time
y_est = f(t)
def fb2(vars,x):
    alpha = vars[0] 
    ceta1 =vars[1] 
    ceta2 =vars[2] 
    sigma1 =vars[3] 
    sigma2 =vars[4]
    
    return (1/pi*(alpha*(arctan(x/ceta1*sigma1) + (1-ceta1)*(sigma1*x)/(x**2+(ceta1*sigma1)**2)) + (1-alpha)*(arctan(x/ceta2*sigma2) + (1-ceta2)*(sigma2*x)/(x**2+(ceta2*sigma2)**2))) + 1/2)

def rho(vars,x):
    alpha = vars[0] 
    ceta1 =vars[1] 
    ceta2 =vars[2] 
    sigma1 =vars[3] 
    sigma2 =vars[4]
    return -1*alpha*((ceta1 - sigma1)/(sigma1**2 + x**2) - 2*x**2*(ceta1 - sigma1)/(sigma1**2 + x**2)**2 + 1/(sigma1*(1 + x**2/sigma1**2)))/pi + (1 - alpha)*((ceta2 - sigma2)/(sigma2**2 + x**2) - 2*x**2*(ceta2 - sigma2)/(sigma2**2 + x**2)**2 + 1/(sigma2*(1 + x**2/sigma2**2)))/pi


def fb(vars,x):
    alpha = vars[0] 
    ceta1 =vars[1] 
    ceta2 =vars[2] 
    sigma1 =vars[3] 
    sigma2 =vars[4]

    return -0.04723951*(11.3130242127366 - 22.6260484254731*x)*exp(-(1.68174196985868 - 3.36348393971735*(1./pi*(alpha*(arctan(x/ceta1*sigma1) + (1-ceta1)*(sigma1*x)/(x**2+(ceta1*sigma1)**2)) + (1-alpha)*(arctan(x/ceta2*sigma2) + (1-ceta2)*(sigma2*x)/(x**2+(ceta2*sigma2)**2))) + 1/2))**2)

#vars =[  0.97732334,   6.18834128,   4.79921026,   6.38894934,  10.15338487]
#vars =[  2.10307809,   4.1307014 ,  12.01390234,  43.14173917,  53.28484069]
vars =[  7.97728922e+00,   4.53711669e-02,   4.01352092e+01,
	         1.82789310e+01,   5.08345263e+02]
x = np.linspace(-10,10,110)
# create plot
plt.plot(x,fb2(vars,x),'o-',label='original data',markersize=5)
#plt.plot(t,y,'.',label='original data',markersize=5)
#plt.plot(t,y_est,'o-',label='estimate',markersize=1)
plt.xlabel('time')
plt.ylabel('sensor readings')
plt.title('least squares fit of degree 6')
plt.savefig('rhosample.png')
