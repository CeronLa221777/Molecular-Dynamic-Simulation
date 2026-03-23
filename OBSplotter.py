import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation
import os

base_name = "obs_3D_NVT_UNI_N200_rho0.050RAD_v1.0_pert_period"
input_path = f"results/{base_name}.dat"
output_path = f"results/{base_name}.png"

if "NVT" in base_name:
    ensemble_label = "NVT Ensemble (Constant Temperature)"
else:
    ensemble_label = "NVE Ensemble (Constant Energy)"

t, K, U, E, T = np.loadtxt(input_path, delimiter=" ", skiprows=1, unpack = True)

#crear la figura
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10, 8), sharex = True)

# Añadir el título general dinámico
fig.suptitle(f'Evolution of the System Observables\n[{ensemble_label}]', fontsize=14, fontweight='bold')

ax1.plot(t, K, label=f'K/N (Kinetic per particle)')
ax1.plot(t, U, label=f'U/N (Potential per particle)')
ax1.plot(t, E, label=f'E/N (Total per particle)')
ax1.set_xlabel('Time (t)')
ax1.set_ylabel('Energy per particle')
ax1.legend()
ax1.grid(True)

ax2.plot(t, T, label='Temperature (T)', color = "purple")
ax2.set_xlabel('Time (t)')
ax2.set_ylabel('Temperature')
ax2.legend()
ax2.grid(True)

plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()




#os.system("cd compound && convert -delay 10 -loop 0 *.png compound.gif && xdg-open compound.gif")