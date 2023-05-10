#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 11:51:43 2023

@author: nico
"""
%reset
import numpy as np
import matplotlib.pyplot as plt
import time
import pylab as p

m = 1 # 9.1e-31
a = 1 # 1e-9
A = np.sqrt(2/a)
hbar = 1 # 6.626e-34

x = np.linspace(0, a, 500)

psi1 = A*np.sin(np.pi/a*x)
psi2 = A*np.sin(2*np.pi/a*x)

E1 = (1**2)*(np.pi**2)*(hbar**2)/(2*m*(a**2))
E2 = (2**2)*(np.pi**2)*(hbar**2)/(2*m*(a**2))

c1 = 1/np.sqrt(2)
c2 = np.sqrt(1 - c1**2)

t = 0

while t <= 1:
    psi = c1*psi1*np.exp(-1j*E1*t/hbar) + c2*psi2*np.exp(-1j*E2*t/hbar)

    plt.plot(np.real(psi))
    plt.plot(np.imag(psi))
    plt.plot(np.abs(psi)**2)
    
    plt.xlim(0,len(psi))
    plt.ylim(-2,3.5)

    plt.draw()
    plt.savefig(str(round(t*100))+'.png')
    plt.pause(0.1)
    plt.clf()
    
    t += 0.01