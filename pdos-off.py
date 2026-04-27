from math import pi, log, exp, sqrt
from scipy.constants import physical_constants
import numpy as np
import matplotlib.pyplot as plt
h_to_ev = physical_constants["Hartree energy in eV"][0]

SIGMA=0.1
STEP=0.1
PADDING = 1.5

data_H = np.loadtxt("ZnS_slab-pdos-k4-1.pdos")
data_Zn = np.loadtxt("ZnS_slab-pdos-k1-1.pdos")

energies = data_H[:, 1]*h_to_ev
coefficients_H = data_H[:, 3:]
coefficients_Zn = data_Zn[:, 3:]
sum_of_coefficients_H = np.sum(coefficients_H, axis=1)
sum_of_coefficients_Zn = np.sum(coefficients_Zn, axis=1)


BINS = int((np.max(energies) - np.min(energies) + 2*PADDING)/STEP)
interval = np.linspace(start=np.min(energies)-PADDING, stop=np.max(energies)+PADDING, num=BINS)


density_H = np.array([])
density_Zn = np.array([])

for b in interval:
    pDOS_H = 0.0
    pDOS_Zn = 0.0
    for e, c2 in zip(energies, sum_of_coefficients_H):
        pDOS_H += (1/(SIGMA*2*pi))*c2*exp((-((b-e)**2)/(2*SIGMA**2)))
    density_H = np.append(density_H, pDOS_H)
    
    for e, c2 in zip(energies, sum_of_coefficients_Zn):
        pDOS_Zn += (1/(SIGMA*2*pi))*c2*exp((-((b-e)**2)/(2*SIGMA**2)))
    density_Zn = np.append(density_Zn, pDOS_Zn)
    

plt.figure(figsize=(8,6))
plt.plot(interval, density_H, color="blue")
plt.plot(interval, density_Zn, color="red")

plt.show()