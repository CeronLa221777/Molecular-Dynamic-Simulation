import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation
import os

t, K, U, E = np.loadtxt("results/obs_3D_UNI_N750_rho0.050RAD_v1.0_pert_period.dat", delimiter=" ", skiprows=1, unpack = True)

plt.figure(figsize=(10, 6))
plt.plot(t, K, label=f'K(t)')
plt.plot(t, U, label=f'U(t)')
plt.plot(t, E, label=f'E(t)')

plt.xlabel('Time (t)')
plt.ylabel('Energy')
plt.title('Evolution of the System Observables')
plt.legend()
plt.grid(True)
plt.savefig(f"observables2DPTchaos.png", dpi=300, bbox_inches='tight')
plt.show()




#os.system("cd compound && convert -delay 10 -loop 0 *.png compound.gif && xdg-open compound.gif")