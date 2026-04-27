#esempio

from math import pi, log, exp, sqrt
from scipy.constants import physical_constants
import numpy as np
import matplotlib.pyplot as plt
from re import compile, search
h_to_ev = physical_constants["Hartree energy in eV"][0]

SIGMA=0.05
STEP=0.05
PADDING = 1.5
PATTERN = compile(r"E\(Fermi\) =\s*(-?\d*\.\d*)")

# apre uno dei file per leggere la Fermi, è indifferente quale
with open("ZnS_slab-pdos-k1-1.pdos", "r") as f:
    for line in f:
        if fermi_found := search(PATTERN, line.strip()):
            fermi = float(fermi_found.group(1))*h_to_ev
            break

data_Zn = np.loadtxt("ZnS_slab-pdos-k1-1.pdos")
data_S = np.loadtxt("ZnS_slab-pdos-k2-1.pdos")
data_O = np.loadtxt("ZnS_slab-pdos-k3-1.pdos")
data_H = np.loadtxt("ZnS_slab-pdos-k4-1.pdos")

energies = data_H[:, 1]*h_to_ev

coefficients_Zn = data_Zn[:, 3:]
coefficients_H = data_H[:, 3:]
coefficients_O = data_O[:, 3:]
coefficients_S = data_S[:, 3:]

sum_of_coefficients_H = np.sum(coefficients_H, axis=1)
sum_of_coefficients_Zn = np.sum(coefficients_Zn, axis=1)
sum_of_coefficients_S = np.sum(coefficients_S, axis=1)
sum_of_coefficients_O = np.sum(coefficients_O, axis=1)

BINS = int((np.max(energies) - np.min(energies) + 2*PADDING)/STEP)
interval = np.linspace(start=np.min(energies)-PADDING, stop=np.max(energies)+PADDING, num=BINS)

density_H = np.array([])
density_Zn = np.array([])
density_O = np.array([])
density_S = np.array([])

for b in interval:
    pDOS_H = 0.0
    pDOS_Zn = 0.0
    pDOS_O, pDOS_S = 0.0, 0.0
    
    for e, c2 in zip(energies, sum_of_coefficients_H):
        pDOS_H += (1/(SIGMA*2*pi))*c2*exp((-((b-e)**2)/(2*SIGMA**2)))
    density_H = np.append(density_H, pDOS_H)
    
    for e, c2 in zip(energies, sum_of_coefficients_Zn):
        pDOS_Zn += (1/(SIGMA*2*pi))*c2*exp((-((b-e)**2)/(2*SIGMA**2)))
    density_Zn = np.append(density_Zn, pDOS_Zn)
    
    for e, c2 in zip(energies, sum_of_coefficients_O):
        pDOS_O += (1/(SIGMA*2*pi))*c2*exp((-((b-e)**2)/(2*SIGMA**2)))
    density_O = np.append(density_O, pDOS_O)
    
    for e, c2 in zip(energies, sum_of_coefficients_S):
        pDOS_S += (1/(SIGMA*2*pi))*c2*exp((-((b-e)**2)/(2*SIGMA**2)))
    density_S = np.append(density_S, pDOS_S)
    

plt.figure(figsize=(8,6))
plt.plot(interval-fermi, density_H, color="blue", label="H")
plt.plot(interval-fermi, density_Zn, color="red", label="Zn")
plt.plot(interval-fermi, density_O, color="cyan", label="O")
plt.plot(interval-fermi, density_S, color="yellow", label="S")
plt.legend()
plt.xlim(-5, 5)
plt.ylim(0,50)
plt.show()
