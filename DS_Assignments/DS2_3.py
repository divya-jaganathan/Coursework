#!/usr/bin/env python

""" 
Dynamical Systems: Assignment #1, Problem 3
Written by: Divya Jagannathan
Date: 13 January, 2019
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib
import matplotlib.pyplot as plt
import math
import sys

# Get rho from call argument

rho=float(sys.argv[1])
print ("rho: %f" %(rho))

# Defining function lorenz

def lorenz(x,t):
 sigma=10.0
 beta=8.0/3.0
 return sigma*(x[1]-x[0]),x[0]*(rho-x[2])-x[1],x[0]*x[1]-beta*x[2];

# Finding fixed points of lorenz

def fixedpts():
 sigma=10.0
 beta=8.0/3.0
 xfplus = (beta*(rho-1.0))**0.5
 zfplus = rho-1.0
 f1=[0,0,0]
 fplus=[xfplus,xfplus,zfplus]
 fminus=[-xfplus,-xfplus,zfplus]
 return f1,fplus,fminus

# Linearized function lorenz

def linlorenz(x,t):
 sigma=10.0
 beta=8.0/3.0
 xfpts=fixedpts()
 xf=xfpts[0]
 return sigma*(x[1]-x[0]), x[0]*(rho-xf[2])-x[1]-xf[0]*x[2],xf[1]*x[0]+xf[0]*x[1]-beta*x[2];

# Initial conditions and time span

xa0=[0.5,0.0,0.0]
xb0=[0.5+10**-10,0.0+10**-10,0.0+10**-10]
t=np.arange(0.0,5.0,0.1)

# Solving ode

fixedpoints=fixedpts()
print(fixedpoints)

x1=odeint(lorenz,xa0,t)
x2=odeint(lorenz,xb0,t)

xl1=odeint(linlorenz,xa0,t)
xl2=odeint(linlorenz,xb0,t)

# L2 Norm
norm=((x1[:,1]-x2[:,1])**2+(x1[:,2]-x2[:,2])**2+(x1[:,0]-x2[:,0])**2)**0.5
normlin=((xl1[:,1]-xl2[:,1])**2+(xl1[:,2]-xl2[:,2])**2+(xl1[:,0]-xl2[:,0])**2)**0.5

# Plotting
plt.subplot(221)
plt.plot(t,xl1[:,0],label='A')
plt.plot(t,xl2[:,0],label='B')
plt.xlabel('$t$')
plt.ylabel('$x_1(t)$')
plt.legend(loc='upper left')

plt.subplot(222)
plt.plot(t,xl1[:,1],label='A')
plt.plot(t,xl2[:,1],label='B')
plt.xlabel('$t$')
plt.ylabel('$x_2(t)$')
plt.legend(loc='upper left')

plt.subplot(223)
plt.plot(t,xl1[:,2],label='A')
plt.plot(t,xl2[:,2],label='B')
plt.xlabel('$t$')
plt.ylabel('$x_3(t)$')

plt.subplot(224)
plt.semilogy(t,normlin)
plt.xlabel('$t$')
plt.ylabel('$||x_A(t)-x_B(t)||_2$')
plt.savefig("DS_3a_350.png")

plt.title('Linearized Lorenz for rho=350')

