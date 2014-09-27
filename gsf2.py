#!/usr/bin/env python
# Date : 2011-01-14 12:23
# Author : zenghui liu  akakcolin(at)163.com
"""
some symbol calculations 
"""
from __future__ import division
import numpy as np
from sympy import *
x = symbols('x')
b = symbols('b')
xi,xi1,xi2 =symbols('xi','xi1','xi2')
Kb,alpha,ceta1,ceta2,sigma1,sigma2 = symbols('Kb','alpha','ceta1','ceta2','sigma1','sigma2')
gaussian = map(Function,'gaussian')
def test(x):
    return atan(x) 
fb= map(Function,'fb')

def fb(x):
    return b/pi*(alpha*(atan((x-xi1)/ceta1*sigma1) + (1-ceta1)*(sigma1*(x-xi1))/((x-xi1)**2+(ceta1*sigma1)**2)) + (1-alpha)*(atan((x-xi2)/ceta2*sigma2) + (1-ceta2)*(sigma2*(x-xi2))/((x-xi2)**2+(ceta2*sigma2)**2))) + b/2

def fb2(x):   
    return 1/pi*((alpha*(atan(x/sigma1) + (ceta1-sigma1)*x)/(x**2+sigma1**2)) + (1-alpha)*(atan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2))) + 1/2

def fb2_1(x):
    return 1/pi*(alpha*(atan(x/sigma1) + (ceta1-sigma1)*x/(x**2+sigma1**2)))

def fb2_2(x):
    return 1/pi*((1-alpha)*(atan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2)))
def gaussians(x):
    return -0.2 +(1.35/(0.43*sqrt(pi/2)))*exp(-2*((x-0.5)/0.43)**2)
    #return -0.0085555 +(0.03441615/(0.50282402*sqrt(pi/2)))*exp(-2*((x-0.5)/0.50282402)**2)

#rho1 = Derivative(fb2_1(xi),xi)
#rho2 = Derivative(fb2_2(xi),xi)
#print rho1.doit()
#print rho2.doit()

def rho1(xi):
    return (alpha*((ceta1 - sigma1)/(sigma1**2 + xi**2) - 2*xi**2*(ceta1 - sigma1)/(sigma1**2 + xi**2)**2 + 1/(sigma1*(1 + xi**2/sigma1**2)))/pi+ (1 - alpha)*((ceta2 - sigma2)/(sigma2**2 + xi**2) - 2*xi**2*(ceta2 - sigma2)/(sigma2**2 + xi**2)**2 + 1/(sigma2*(1 + xi**2/sigma2**2)))/pi)/(x-xi)

def rhotest(x):
    return 1/2*Kb*sin((1/pi*((alpha*(atan(x/sigma1) + (ceta1-sigma1)*x)/(x**2+sigma1**2)) + (1-alpha)*(atan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2))) + 1/2)*4*pi)
#rho1out=series(rhotest(x),x,0,3)
#print rho1out.doit()

def rhoseries(xi):
    return ceta2/(pi*sigma2**2*x) + alpha*ceta1/(pi*sigma1**2*x) - alpha*ceta2/(pi*sigma2**2*x) + ceta2*xi/(pi*sigma2**2*x**2) + alpha*ceta1*xi/(pi*sigma1**2*x**2) - alpha*ceta2*xi/(pi*sigma2**2*x**2) + 2*xi**2/(pi*sigma2**3*x) + ceta2*xi**2/(pi*sigma2**2*x**3) - 3*ceta2*xi**2/(pi*sigma2**4*x) - 2*alpha*xi**2/(pi*sigma2**3*x) + 2*alpha*xi**2/(pi*sigma1**3*x) + alpha*ceta1*xi**2/(pi*sigma1**2*x**3) - alpha*ceta2*xi**2/(pi*sigma2**2*x**3) - 3*alpha*ceta1*xi**2/(pi*sigma1**4*x) + 3*alpha*ceta2*xi**2/(pi*sigma2**4*x) + 2*xi**3/(pi*sigma2**3*x**2) + ceta2*xi**3/(pi*sigma2**2*x**4) - 3*ceta2*xi**3/(pi*sigma2**4*x**2) - 2*alpha*xi**3/(pi*sigma2**3*x**2) + 2*alpha*xi**3/(pi*sigma1**3*x**2) + alpha*ceta1*xi**3/(pi*sigma1**2*x**4) - alpha*ceta2*xi**3/(pi*sigma2**2*x**4) - 3*alpha*ceta1*xi**3/(pi*sigma1**4*x**2) + 3*alpha*ceta2*xi**3/(pi*sigma2**4*x**2) 

#rho3=series(rho1(xi),xi,0,4)
#print rho3.doit()
rhointegrate=integrate(rhoseries(xi),(xi,0,1))
print rhointegrate.doit()

def right(x):
    return 1/2*Kb*sin((1/pi*((alpha*(atan(x/sigma1) + (ceta1-sigma1)*x)/(x**2+sigma1**2)) + (1-alpha)*(atan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2))) + 1/2)*4*pi)-(-0.0085555 +(0.03441615/(0.50282402*sqrt(pi/2)))*exp((-2*((1/pi*((alpha*(atan(x/sigma1) + (ceta1-sigma1)*x)/(x**2+sigma1**2)) + (1-alpha)*(atan(x/sigma2) + (ceta2-sigma2)*x/(x**2+sigma2**2))) + 1/2)-0.5)/0.50282402)**2))
