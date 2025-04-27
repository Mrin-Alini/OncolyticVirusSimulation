import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import Functions_OncolyticVirus as FOV
import Parameters_OncolyticVirus as POV

# Original parameters and settings
y0_base = POV.y0_Ad1  # Base initial conditions with variable V0

# Define a range for initial V0 values
V_range = np.linspace(0.001, 0.02, 20)

# Find the optimal V0 for the earliest tumor volume peak
earliest_V0, min_peak_time, peak_y = FOV.find_earliest_peak(y0_base, POV.params_Ad1, FOV.ideal_t_param, V_range)

y0 = y0_base.copy()
y0[3] = earliest_V0  # Update the initial viral load
earliest_y = odeint(FOV.sim, y0, FOV.t, args=(POV.params_Ad1,))

# Plot the optimal solution
if earliest_V0 == None:
    print(f"The earliest tumor volume peak occurs after t = {min_peak_time}")
else:

    f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)

    LineT, = ax1.plot(FOV.t, earliest_y[:, 0], color="blue", label="LineT")
    LineV, = ax2.plot(FOV.t, earliest_y[:, 3], color="red", label="LineV")

    ax1.set_ylabel('Tumor Volume')
    ax2.set_ylabel('Virus Volume')
    ax2.set_xlabel('Time')
    plt.legend(handles=[LineT, LineV])
    plt.show()
    print(f"Optimal V0: {earliest_V0:.4f} gives the earliest peak time at {min_peak_time:.2f} days"
      f" with a volume of {peak_y:.4f}.")