import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import Functions_OncolyticVirus as FOV
import Parameters_OncolyticVirus as POV

# Solve the system of ODEs
y = odeint(FOV.sim, POV.y0_Ad5, FOV.t, args=(POV.params_Ad5,))

# Plotting
f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True, sharey=False, figsize=(6,8))

LineT, = ax1.plot(FOV.t, y[:,0], color="blue", label="LineT")
LineE, = ax2.plot(FOV.t, y[:,1], color="yellow", label="LineE")
LineI, = ax3.plot(FOV.t, y[:,2], color="green", label="LineI")
LineV, = ax4.plot(FOV.t, y[:,3], color="red", label="LineV")
LineF, = ax5.plot(FOV.t, y[:,4], color="black", label="LineF")

plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

ax1.set_title("Volume of Tumor Treated with Ad5.P19 Virus Over Time")

ax2.set_title("Volume of Cells in Eclipse Phase Over Time")

ax3.set_title("Volume of Infected Cancer Cells Over Time")

ax4.set_title("Volume of Ad5.P19 Virus Over Time")

ax5.set_title("Volume of Interferon Over Time")
ax5.set_xlabel('Time')

#ax1.legend(bbox_to_anchor=(1, 1), handles=[LineT, LineE, LineI, LineV, LineF])

plt.show()

tumor_volume = y[:, 0]
peak_vol_index = np.argmax(tumor_volume)
print (peak_vol_index)
peak_time = FOV.t[peak_vol_index]
peak_vol = tumor_volume[peak_vol_index]
print(peak_time, peak_vol)

print(max(y[:, 0]))
