#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 17:54:42 2024

@author: nico
"""

import numpy as np
import matplotlib.pyplot as plt

# Define k from 2 to 200 with step size 2 (even numbers)
# odd values of k result in coeff = 0
k = np.arange(2, 201, 2)

# Calculate H'kn values
Hpkn = np.sin((k - 1) / 2 * np.pi) / (k - 1) - np.sin((k + 1) / 2 * np.pi) / (k + 1)

# Plot H'kn values
plt.figure()
plt.plot(k, Hpkn, '-*')
plt.xlim([0, 100])
plt.xlabel("k")
plt.ylabel("<k⁰|H'|1⁰>")
plt.show()

# Compute the sum of |Hpkn|² / (1 - k²) up to each value of k
suma = np.zeros(len(Hpkn))
for i in range(len(Hpkn)):
    suma[i] = np.sum(Hpkn[:i+1]**2 / (1 - k[:i+1]**2))

# Evolution of the sum
plt.figure()
plt.plot(k, suma)
plt.xlim([0, 50])
plt.xlabel("k")
plt.ylabel("suma")
plt.show()