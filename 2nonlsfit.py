#!/usr/bin/env python
# Date : 2011-01-14 14:30
# Author : zenghui liu  akakcolin(at)163.com

import matplotlib.pyplot as plt
from numpy import exp,sin,linspace,pi,arctan,sqrt
from lmfit import minimize,Parameters

def fb2(params,x):
    alpha = params['alpha'].value
    ceta1 = params['ceta1'].value
    ceta2 = params['ceta2'].value
    sigma1 =params['sigma2'].value
    sigma2 =params['sigma2'].value

    return (1/pi*(alpha*(arctan(x/ceta1*sigma1) + (1-ceta1)*(sigma1*x)/(x**2+(ceta1*sigma1)**2)) + (1-alpha)*(arctan(x/ceta2*sigma2) + (1-ceta2)*(sigma2*x)/(x**2+(ceta2*sigma2)**2))) + 1./2)

def gaussian(params,x):
    alpha = params['alpha'].value
    ceta1 = params['ceta1'].value 
    ceta2 = params['ceta2'].value 
    sigma1 =params['sigma2'].value 
    sigma2 =params['sigma2'].value 

    return -0.1371 +(0.5514/(0.5028*sqrt(pi/2)))*exp((-2*((1/pi*((alpha*(arctan(x/sigma1) + (ceta1-sigma1)*x)/(x**2+sigma1**2)) + (1-alpha)*(arctan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2))))/0.5028)**2))

def rho(params,x):
    alpha = params['alpha'].value
    ceta1 = params['ceta1'].value
    ceta2 = params['ceta2'].value
    sigma1 =params['sigma2'].value
    sigma2 =params['sigma2'].value

    return 1./pi*(alpha*sigma1*(1 - ceta1)/(x**2 + ceta1**2*sigma1**2) - 2*sigma1*x**2*(1 - ceta1)/(x**2 + ceta1**2*sigma1**2)**2 + sigma1/(ceta1*(1 + sigma1**2*x**2/ceta1**2))+(1 - alpha)*(sigma2*(1 - ceta2)/(x**2 + ceta2**2*sigma2**2) - 2*sigma2*x**2*(1 - ceta2)/(x**2 + ceta2**2*sigma2**2)**2 + sigma2/(ceta2*(1 + sigma2**2*x**2/ceta2**2)))) 

def objfun(params,x,data): 
    alpha = params['alpha'].value
    ceta1 = params['ceta1'].value
    ceta2 = params['ceta2'].value
    sigma1 =params['sigma2'].value
    sigma2 =params['sigma2'].value
    Kb = params['Kb'].value

   # model = Kb/2/pi*sin((1/pi*((alpha*(arctan(x/sigma1) + (ceta1-sigma1)*x)/(x**2+sigma1**2)) + (1-alpha)*(arctan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2))) + 1/2)*4*pi)
  #  return (model-gaussian(params,x))

    return  Kb/2/pi*sin((1/pi*((alpha*(arctan(x/sigma1) + (ceta1-sigma1)*x)/(x**2+sigma1**2)) + (1-alpha)*(arctan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2))) + 1/2)*4*pi)+0.1371 -(0.5514/(0.5028*sqrt(pi/2)))*exp((-2*((1/pi*((alpha*(arctan(x/sigma1) + (ceta1-sigma1)*x)/(x**2+sigma1**2)) + (1-alpha)*(arctan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2))))/0.5028)**2))
    


params = Parameters()
params.add('alpha',value=0.98984467)
params.add('ceta1',value=76.79684648)
params.add('ceta2',value=30.77261877)
params.add('sigma1',value=77.1949924)
params.add('sigma2',value=70.06017928)
params.add('Kb',value=24.74401521,vary=False)
for name,par in params.items():
    print ' %s =%.4f +/- ' %(name, par.value)

x =linspace(-10,10,101)
data =linspace(0,0,101)
out = minimize(objfun,params,args=(x,data),engine='leastsq')
print out.redchi
print "Best-Fit Values:"
for name,par in params.items():
    print ' %s =%.4f +/- %.4f' %(name, par.value,par.stderr)
print params

#plt.plot(x,gau(x),'-',label='dgoriginal data',markersize=1)
#plt.plot(x,fb2(params,x),'o-',label='EAM',markersize=5)
plt.plot(x,rho(params,x),'-',label='EAM',markersize=5)
#plt.plot(t,y,'.',label='original data',markersize=5)
#plt.plot(t,y_est,'o-',label='estimate',markersize=1)
plt.xlabel('x/b',fontsize=15)
plt.ylabel('rho(x)/b',fontsize=15)
#plt.ylabel('fb(x)/b',fontsize=15)
plt.title('dislocation density rho(x)',fontsize=15)
#plt.title('dislocation disregistry fp(x)',fontsize=15)
plt.legend(loc='best')
plt.fill_between(x,0,rho(params,x),facecolor='green',alpha=0.2)
plt.savefig('rhosample.png')

