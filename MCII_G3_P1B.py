import numpy as np
import matplotlib.pyplot as plt

# tomando n = 1
ncoef = 20
k = np.arange(2, ncoef, 2)  # los k impares dan coef = 0
Hpkn = np.sin((k-1)/2*np.pi)/(k-1) - np.sin((k+1)/2*np.pi)/(k+1)  # elementos H'

# se muestran los sucesivos valores de H'_k1
plt.figure()
plt.plot(k, Hpkn, '-*')
plt.xlim([0, ncoef])
plt.xlabel('k')
plt.ylabel('<k⁰|Hp|1⁰>')
plt.show()

suma = []
# se calcula la suma de |H'k1|²/(1-k²) hasta cada valor de k
for i in range(len(Hpkn)):
    suma.append(np.sum(Hpkn[:i+1]**2 / (1 - k[:i+1]**2)))
    print(f"coef: {Hpkn[i]**2 / (1 - k[i]**2)} - suma: {suma[i]}")

# se muestra la evolución de la suma
plt.figure()
plt.plot(k, suma)
plt.xlabel('k')
plt.ylabel('suma')
plt.show()

# tomando n = 2
k = np.arange(1, ncoef, 2)  # los k pares dan coef = 0
Hpkn = np.sin((k-2)/2*np.pi)/(k-2) - np.sin((k+2)/2*np.pi)/(k+2)  # elementos H'

# se muestran los sucesivos valores de H'_k1
plt.figure()
plt.plot(k, Hpkn, '-*')
plt.xlabel('k')
plt.ylabel('<k⁰|Hp|2⁰>')
plt.show()

suma = []
# se calcula la suma de |H'k2|²/(4-k²) hasta cada valor de k
for i in range(len(Hpkn)):
    suma.append(np.sum(Hpkn[:i+1]**2 / (4 - k[:i+1]**2)))
    print(f"coef: {Hpkn[i]**2 / (4 - k[i]**2)} - suma: {suma[i]}")

# se muestra la evolución de la suma
plt.figure()
plt.plot(k, suma)
plt.xlabel('k')
plt.ylabel('suma')
plt.show()
