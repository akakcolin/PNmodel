#!/usr/bin/env python
# Date : 2011-01-14 15:23
# Author : zenghui liu  akakcolin(at)163.com
"""
This file is used to fit the GSF energy curve to a gaussian form function

The unit must be considered right.
"""
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from lmfit import minimize,Parameters

# read data from gsf data  (format: x y, each line)
def read_data(filename):
    """
    read data from file
    """
    f = open(filename,'r')
    xlist=[]
    ylist=[]
    allpair=[]
    for line in f:
	pair = line.split()
	if pair == "None":
	    pass
	else: 
	    x = float(pair[0]) 
	    y = float(pair[1]) 
	    allpair.append([x,y])
    lastpair = np.array(allpair)
    return lastpair 

dat = read_data('gsf.dat')

# define the bugers 
bugers=2.863946

# the data type shoud be float ?
x = np.array([i/bugers for i in dat[:,0]])
data = np.array([i for i in dat[:,1]])

# GSF is fitted to a gaussian function
# the objective function to fit
def residual(params,x,data):
    y0 = params['y0'].value
    mu = params['mu'].value
    sigma = params['sigma'].value
    A = params['A'].value
    
    model = y0 +(A/(sigma*np.sqrt(np.pi/2)))*np.exp(-2*((x-mu)/sigma)**2) 

    return (data-model)

def gaussians(params,x):
    y0 = params['y0'].value
    mu = params['mu'].value
    sigma = params['sigma'].value
    A = params['A'].value
    
    model = y0 +(A/(sigma*np.sqrt(np.pi/2)))*np.exp(-2*((x-mu)/sigma)**2)

    return model

params = Parameters()
params.add('y0',value=0)
params.add('mu',value=0.5,min=0.0)
params.add('sigma',value=0.2)
params.add('A',value=3.0)

vars=[0.5,0.50282405,0.03441615, 0.000]

#[outs,rest] = leastsq(residual,vars,args=(x,data))
#print outs
out = minimize(residual,params,args=(x,data))
print out.chisqr
print "Best-Fit Values:"
for name,par in params.items():
        print ' %s =%.4f +/- %.4f' %(name, par.value,par.stderr)

#0.00296921993569
#Best-Fit Values: 
#    y0 =-0.1371 +/- 0.0264 
#    mu =0.5000 +/- 0.0025 
#    sigma =0.5028 +/- 0.0167 
#    A =0.5514 +/- 0.0321

xtem = np.linspace(-0.0,1.0,100)
y_est=gaussians(params,xtem)

#create plot
#plt.plot(xtem,y_deriv(xtem),'o-',label='estimate',markersize=1)
plt.plot(xtem,y_est,'-',label='estimate',markersize=1)
plt.plot(x,data,'.',label='original',markersize=5)
plt.xlabel('x/burger',fontsize=15)
plt.ylabel(' GSF energy for [1-11](110) J/m2',fontsize=15)
plt.title('least squares fit for GSF',fontsize=15)
plt.savefig('gsffit2.png')
