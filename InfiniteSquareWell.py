import numpy as np
import matplotlib.pyplot as plt
import time
import pylab as p

# constants
m = 1 # 9.1e-31
a = 1 # 1e-9
A = np.sqrt(2/a)
hbar = 1 # 6.626e-34

# positions
x = np.linspace(0, a, 500)

# stationary states
psi1 = A*np.sin(np.pi/a*x)
psi2 = A*np.sin(2*np.pi/a*x)
# add as many as desired, with psin = A*np.sin(n*np.pi/a*x)

# energy eigenvalues of stationary states
E1 = (1**2)*(np.pi**2)*(hbar**2)/(2*m*(a**2))
E2 = (2**2)*(np.pi**2)*(hbar**2)/(2*m*(a**2))
# add as many as desired, with En = (n**2)*(np.pi**2)*(hbar**2)/(2*m*(a**2))

# coefficients of the expansion Psi(x,0) = sum_i c_i psi_i for the initial state Psi(x,0)
c1 = 1/np.sqrt(2)
c2 = np.sqrt(1 - c1**2)

# time evolution
tsim = 10 # simulation time 
dt = 0.01 # timestep
while t <= tsim:
    # Psi(x,t)
    psi = c1*psi1*np.exp(-1j*E1*t/hbar) + c2*psi2*np.exp(-1j*E2*t/hbar)

    plt.plot(np.real(psi))    # plot real part of wavefunction
    plt.plot(np.imag(psi))    # plot imaginary part of wavefunction
    plt.plot(np.abs(psi)**2)  # plot square modulus of wavefunction
    
    plt.xlim(0,len(psi)) # set x limits
    plt.ylim(-2,3.5)     # set y limites

    # draw, (save if wanted), pause, clear
    plt.draw()
    # plt.savefig(str(round(t*100))+'.png')
    plt.pause(0.1)
    plt.clf()
    
    # advance time
    t += dt
