#!/usr/bin/env python

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

# Initial conditions and time span

xa0=[1.0,1.0,1.0]
xb0=[1.0+10**-10,1.0+10**-10,1.0+10**-10]
t=np.arange(0.0,100.0,0.1)

# Solving ode

x1=odeint(lorenz,xa0,t)
x2=odeint(lorenz,xb0,t)

# L2 Norm
norm=((x1[:,1]-x2[:,1])**2+(x1[:,2]-x2[:,2])**2+(x1[:,0]-x2[:,0])**2)**0.5

# Plotting

plt.figure(figsize=(11,11))

plt.subplot(221)
plt.plot(t,x1[:,0],label='A')
plt.plot(t,x2[:,0],label='B')
plt.xlabel('$t$')
plt.ylabel('$x_1(t)$')
plt.legend(loc='upper left')

plt.subplot(222)
plt.plot(t,x1[:,1],label='A')
plt.plot(t,x2[:,1],label='B')
plt.xlabel('$t$')
plt.ylabel('$x_2(t)$')
plt.legend(loc='upper left')

plt.subplot(223)
plt.plot(t,x1[:,2],label='A')
plt.plot(t,x2[:,2],label='B')
plt.xlabel('$t$')
plt.ylabel('$x_3(t)$')

plt.subplot(224)
plt.semilogy(t,norm)
plt.xlabel('$t$')
plt.ylabel('$||x_A(t)-x_B(t)||_2$')
plt.savefig("DS_1b.png")

