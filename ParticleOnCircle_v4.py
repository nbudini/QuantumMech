#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 27 13:41:42 2025

@author: nico
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import sys

# Define constants (setting them to 1 to capture the physics)
m = 1        # mass
hbar = 1     # Planck's constant
A = 1/np.sqrt(2*np.pi)  # normalization constant

# Define range for polar angle and grid size
npts = 1000
theta = np.linspace(0, 2*np.pi, npts)

# Define eigenstates as a function of quantum number n and polar angle
def psi_n_theta(n, theta):
    return A * np.exp(1j * n * theta)

# Define initial wavefunction Psi(theta,0) as a Gaussian
sigma = np.pi / 200
xcen = np.pi

def Psi0(x):
    return 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-(x - xcen)**2 / (2*sigma**2))

# Calculate normalization factor
normfactor, _ = quad(lambda x: np.abs(Psi0(x))**2, 0, 2*np.pi)
def Psi0n(x):
    return Psi0(x) / np.sqrt(normfactor)

# Get normalized coefficients of the eigenstates expansion
ncoeff = 20

def get_coeffs(ncoeff, Psi0_func, psi_n_theta_func):
    coeffs = []
    for n in range(1, ncoeff+1):
        integrand = lambda theta: np.conj(psi_n_theta_func(n, theta)) * Psi0_func(theta)
        real_part = quad(lambda theta: np.real(integrand(theta)), 0, 2*np.pi)[0]
        imag_part = quad(lambda theta: np.imag(integrand(theta)), 0, 2*np.pi)[0]
        coeff = real_part + 1j * imag_part
        coeffs.append(coeff)
    return np.array(coeffs)

coeff = get_coeffs(ncoeff, Psi0n, psi_n_theta)

# Calculate expanded wavefunction at t=0
Psi0_exp = np.zeros_like(theta, dtype=complex)
for i, c in enumerate(coeff, start=1):
    Psi0_exp += c * psi_n_theta(i, theta)
Psi0_exp /= np.sqrt(np.sum(np.abs(Psi0_exp)**2))

# Plot initialization
fig = plt.figure(figsize=(10, 8))
plt.subplots_adjust(hspace=0.5)

ax1 = plt.subplot(221)
h1, = ax1.plot(theta, np.real(Psi0_exp), linewidth=2)
ax1.set_xlim([0, 2*np.pi])
ax1.set_ylim([-0.15, 0.17])
ax1.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
ax1.set_xticklabels(['0', 'π/2', 'π', '3π/2', '2π'])
ax1.set_ylabel('Re(Ψ)')

ax2 = plt.subplot(222)
h2, = ax2.plot(theta, np.imag(Psi0_exp), 'r', linewidth=2)
ax2.set_xlim([0, 2*np.pi])
ax2.set_ylim([-0.15, 0.17])
ax2.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
ax2.set_xticklabels(['0', 'π/2', 'π', '3π/2', '2π'])
ax2.set_ylabel('Im(Ψ)')

ax3 = plt.subplot(223)
h3, = ax3.plot(theta, np.abs(Psi0_exp)**2, 'k', linewidth=2)
ax3.set_xlim([0, 2*np.pi])
ax3.set_ylim([0, 0.025])
ax3.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
ax3.set_xticklabels(['0', 'π/2', 'π', '3π/2', '2π'])
ax3.set_xlabel('θ coordinate')
ax3.set_ylabel('|Ψ|^2')

# Polar plot
ax4 = plt.subplot(224, polar=True)
r = np.ones_like(theta)
c = np.abs(Psi0_exp)**2
h4 = ax4.scatter(theta, r, c=c, s=50, cmap='turbo')
ax4.set_rlim(0, 1.1)
ax4.set_yticklabels([])
ax4.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
ax4.set_xticklabels(['0', 'π/2', 'π', '3π/2'])

# Colorbar
plt.colorbar(h4, ax=ax4, orientation='vertical')

# Time evolution
tf = 12.565
dt = 0.005
times = np.arange(0, tf, dt)
r_val = 1  # constant radial value

for t in times:
    Psi_t = np.zeros_like(theta, dtype=complex)
    for j, c in enumerate(coeff, start=1):
        E_n = (j**2 * hbar**2) / (2 * m * r_val**2)
        Psi_t += c * psi_n_theta(j, theta) * np.exp(-1j * E_n * t / hbar)
    Psi_t /= np.sqrt(np.sum(np.abs(Psi_t)**2))

    # Update plots
    h1.set_ydata(np.real(Psi_t))
    h2.set_ydata(np.imag(Psi_t))
    h3.set_ydata(np.abs(Psi_t)**2)
    h4.set_array(np.abs(Psi_t)**2)

    plt.suptitle(f't = {t:.3f}', fontsize=14)
    plt.pause(0.01)

    # Uncomment to save frames as png images
    
    # filename = f'frame_{int(t/dt):04d}.png'
    # plt.savefig(filename, dpi=300)
    
    # To create mp4 video from these png images use this command afterwards
    # in a Linux terminal:
    # $ ffmpeg -framerate 20 -i frame_%04d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p output.mp4
    
    # (Install ffmpeg with: $ sudo apt install ffmpeg)

plt.show()
