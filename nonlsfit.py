#!/usr/bin/env python
# Date : 2011-01-14 11:32
# Author : zenghui liu  akakcolin(at)163.com
import numpy as np
from numpy import exp,sin,pi,arctan,linspace,sqrt
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def objfun(vars,x): 
    
    alpha = vars[0]
    ceta1 =vars[1]
    ceta2 =vars[2]
    sigma1 =vars[3]
    sigma2 =vars[4]
    Kb = vars[5]
    


    return  Kb/2/pi*sin((1/pi*((alpha*(arctan(x/sigma1) + (ceta1-sigma1)*x)/(x**2+sigma1**2)) + (1-alpha)*(arctan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2))) + 1/2)*4*pi)+0.02 -(1.35/(0.43*sqrt(pi/2)))*exp((-2*((1/pi*((alpha*(arctan(x/sigma1) + (ceta1-sigma1)*x)/(x**2+sigma1**2)) + (1-alpha)*(arctan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2))))/0.43)**2))

#vars=[  0.97719638,   6.18822729,   4.86484639,   6.38904662,  10.22759488]
#vars=[  0.98982043,  76.79144964 , 31.27609979 , 77.20035487 , 69.5375197,20 ]
vars= [0.98984467 , 76.79684648 , 30.77261877 , 77.1949924  , 70.06017928, 24.74401521]
#vars=[ 2.,   6.,   12.6,  47.5285876 ,  52.47275106]
x = linspace(-5,5,80)
[out,req]= leastsq(objfun,vars,args=x)
print out

def fb2(vars,x):
    alpha = vars[0]
    ceta1 =vars[1]
    ceta2 =vars[2]
    sigma1 =vars[3]
    sigma2 =vars[4]

    return (1/pi*(alpha*(arctan(x/ceta1*sigma1) + (1-ceta1)*(sigma1*x)/(x**2+(ceta1*sigma1)**2)) + (1-alpha)*(arctan(x/ceta2*sigma2) + (1-ceta2)*(sigma2*x)/(x**2+(ceta2*sigma2)**2))) + 1./2)

def rho(vars,x):
    alpha = vars[0]
    ceta1 =vars[1]
    ceta2 =vars[2]
    sigma1 =vars[3]
    sigma2 =vars[4]
    return 1./pi*(alpha*sigma1*(1 - ceta1)/(x**2 + ceta1**2*sigma1**2) - 2*sigma1*x**2*(1 - ceta1)/(x**2 + ceta1**2*sigma1**2)**2 + sigma1/(ceta1*(1 + sigma1**2*x**2/ceta1**2))+(1 - alpha)*(sigma2*(1 - ceta2)/(x**2 + ceta2**2*sigma2**2) - 2*sigma2*x**2*(1 - ceta2)/(x**2 + ceta2**2*sigma2**2)**2 + sigma2/(ceta2*(1 + sigma2**2*x**2/ceta2**2))) ) 

def gaussian(x):
    return -0.2 +(1.35/(0.43*sqrt(pi/2)))*exp(-2*((x-0.5)/0.43)**2)
def derivative(vars,f):
    """
    compute the numerical derivative 
    """
    def df(vars,x,h=0.1e-5):
	return (f(vars,x+h/2)-f(vars,x-h/2))/h
    return df
def derivative2(f):
    """
    compute the numerical derivative 
    """
    def df(x,h=0.1e-5):
	return (f(x+h/2)-f(x-h/2))/h
    return df

def fb(vars,x):
    alpha = vars[0]
    ceta1 =vars[1]
    ceta2 =vars[2]
    sigma1 =vars[3]
    sigma2 =vars[4]

    return -0.04723951*(11.3130242127366 - 22.6260484254731*x)*exp(-(1.68174196985868 - 3.36348393971735*(1./pi*(alpha*(arctan(x/ceta1*sigma1) + (1-ceta1)*(sigma1*x)/(x**2+(ceta1*sigma1)**2)) + (1-alpha)*(arctan(x/ceta2*sigma2) + (1-ceta2)*(sigma2*x)/(x**2+(ceta2*sigma2)**2))) + 1/2))**2)

#vars =[  0.97732334,   6.18834128,   4.79921026,   6.38894934,  10.15338487]
#vars =[  2.10307809,   4.1307014 ,  12.01390234,  43.14173917,  53.28484069]

#vars =[  7.97728922e+00,   4.53711669e-02,   4.01352092e+01,
x = np.linspace(-10,10,110)
# create plot
gau = derivative2(gaussian)
dg = derivative(out,fb2)
#print dg(out,0)
print fb2(vars,1)
plt.plot(x,gau(x),'-',label='dgoriginal data',markersize=1)
#plt.plot(x,fb2(out,x),'o-',label='original data',markersize=5)
#plt.plot(x,rho(out,x),'o-',label='original data',markersize=5)
#plt.plot(t,y,'.',label='original data',markersize=5)
#plt.plot(t,y_est,'o-',label='estimate',markersize=1)
plt.xlabel('x/b')
plt.ylabel('rho')
plt.title('least squares fit of degree 6')
plt.savefig('rhosample.png')
