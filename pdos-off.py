from math import pi, log, exp, sqrt
from scipy.constants import physical_constants
import numpy as np
import matplotlib.pyplot as plt
from re import compile, search

plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 15,
    "ytick.major.size": 0,
    "ytick.labelsize": 0,
    })
    
COLORS = {
    "Zn": "#141410",
    "O": "#e64727",
    "H": "#ebe1df",
    "Cu": "#4e86d9",
    "S": "#dbdb14",
    "Se": "#eb8f34"
    }

# constants definitions
H_TO_EV = physical_constants["Hartree energy in eV"][0]
SIGMA=0.1
STEP=0.03
PADDING = 1.5
FERMI_PATTERN = compile(r"E\(Fermi\) =\s*(-?\d*\.\d*)")

# path definitions
PATH1="C:\\Users\\franc\\Desktop\\test-python\\DOS\\DOMENICO\\ZnO\\dopato\\OH-\\"
#PATH2="C:\\Users\\franc\\Desktop\\test-python\\DOS\\DOMENICO\\ZnSe\\dopato\\H2O+"

# functions definitions
def calculate_DOS(folder_path, file_name, require_pdos=False, fermi_correction=False):
    
    if not fermi_correction:
        fermi = 0.0
    else:
        with open(f"{folder_path}{file_name}", "r") as f:
            for line in f:
                if fermi_found := search(FERMI_PATTERN, line.strip()):
                    fermi = float(fermi_found.group(1))*H_TO_EV
                    break
            
    data = np.loadtxt(f"{folder_path}{file_name}")
    energies = data[:, 1]*H_TO_EV
    if require_pdos:
        coefficients = data[:, 3:]
        sum_of_coefficients = np.sum(coefficients, axis=1)
    else:
        coefficients = sum_of_coefficients = np.ones(len(energies))
    
    BINS = int((np.max(energies) - np.min(energies) + 2*PADDING)/STEP)
    interval = np.linspace(start=np.min(energies)-PADDING, stop=np.max(energies)+PADDING, num=BINS)
    density = np.array([])
    for b in interval:
        DOS = 0.0
        for e, c2 in zip(energies, sum_of_coefficients):
            DOS += (1/(SIGMA*2*pi))*c2*exp((-((b-e)**2)/(2*SIGMA**2)))
        density = np.append(density, DOS)
        
    #np.savetxt(f"{folder_path}{file_name.split(".")[0]}.dat", np.stack((interval, density), axis=1))

    return interval-fermi, density, fermi
    
i_Zna, d_Zna, f_Zna = calculate_DOS(folder_path=PATH1, file_name="Zna.pdos", require_pdos=True, fermi_correction=True)
i_Znb, d_Znb, f_Znb = calculate_DOS(folder_path=PATH1, file_name="Znb.pdos", require_pdos=True, fermi_correction=True)
i_Cua, d_Cua, f_Cua = calculate_DOS(folder_path=PATH1, file_name="Cua.pdos", require_pdos=True, fermi_correction=True)
i_Cub, d_Cub, f_Cub = calculate_DOS(folder_path=PATH1, file_name="Cub.pdos", require_pdos=True, fermi_correction=True)
i_Oa, d_Oa, f_Oa = calculate_DOS(folder_path=PATH1, file_name="Oa.pdos", require_pdos=True, fermi_correction=True)
i_Ob, d_Ob, f_Ob = calculate_DOS(folder_path=PATH1, file_name="Ob.pdos", require_pdos=True, fermi_correction=True)
i_Ha, d_Ha, f_Ha = calculate_DOS(folder_path=PATH1, file_name="Ha.pdos", require_pdos=True, fermi_correction=True)
i_Hb, d_Hb, f_Hb = calculate_DOS(folder_path=PATH1, file_name="Hb.pdos", require_pdos=True, fermi_correction=True)
# i_Sa, d_Sa, f_Sa = calculate_DOS(folder_path=PATH1, file_name="Sa.pdos", require_pdos=True, fermi_correction=True)
# i_Sb, d_Sb, f_Sb = calculate_DOS(folder_path=PATH1, file_name="Sb.pdos", require_pdos=True, fermi_correction=True)
# i_Sea, d_Sea, f_Sea = calculate_DOS(folder_path=PATH1, file_name="Sea.pdos", require_pdos=True, fermi_correction=True)
# i_Seb, d_Seb, f_Seb = calculate_DOS(folder_path=PATH1, file_name="Seb.pdos", require_pdos=True, fermi_correction=True)

plt.figure(figsize=(8,8))
#pop alpha
plt.plot(i_Ha, d_Ha, color=COLORS["H"], label="H")
# plt.plot(i_Zna, d_Zna, color=COLORS["Zn"], label="Zn")
# plt.plot(i_Cua, d_Cua, color=COLORS["Cu"], label="Cu")
# plt.plot(i_Oa, d_Oa, color=COLORS["O"], label="O")
# plt.plot(i_Sa, d_Sa, color=COLORS["S"], label="S")
# plt.plot(i_Sea, d_Sea, color=COLORS["Se"], label="Se")

#pop beta
plt.plot(i_Hb, -d_Hb, color=COLORS["H"], label="H")
# plt.plot(i_Znb, -d_Znb, color=COLORS["Zn"], label="Zn")
# plt.plot(i_Cub, -d_Cub, color=COLORS["Cu"], label="Cu")
# plt.plot(i_Ob, -d_Ob, color=COLORS["O"], label="O")
# plt.plot(i_Sb, -d_Sb, color=COLORS["S"], label="S")
# plt.plot(i_Seb, -d_Seb, color=COLORS["Se"], label="Se")

# plt.xlim(-5, 5)
# plt.ylim(-32, 32)   

plt.xlabel("energy (eV)")
plt.ylabel("pDOS (a.u.)")
plt.legend()
# plt.axhline(0, linewidth=0.5, color="black")
# plt.savefig(f"{PATH1}pDOS.svg")
plt.show()
