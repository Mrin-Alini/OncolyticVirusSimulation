import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import Functions_OncolyticVirus as FOV
import Parameters_OncolyticVirus as POV

# Original parameters and settings
y0_base = POV.y0_Ad6  # Base initial conditions with variable V0

# Define a range for V0 and e values
V_range = np.linspace(0.001, 0.02, 20)
e_range = np.linspace(1, 60000, 20)
c_range = np.linspace(0.04, 0.3, 20)
alpha_range = np.linspace(4.854*10**-13, 2.421*10**-7, num=20)

# Find the optimal V0 for the earliest tumor volume peak
optimum_params, optimum_peak, acceptable_params_range, acceptable_params, min_peak_time \
    = FOV.find_optimum_and_range_of_many_param(y0_base, POV.params_Ad6,
                                               FOV.ideal_t_param, FOV.t, V_range, c_range, alpha_range)

#optimum_params, optimum_peak, acceptable_params_range, acceptable_params, \
#min_peak_time = FOV.find_optimum_and_range_of_multiple_param(y0_base, POV.params_Ad1,
                                                            #FOV.ideal_t_param, FOV.t, e_range, c_range, alpha_range)
optimum_V0 = optimum_params[0]
#optimum_e = optimum_params[1]
optimum_c = optimum_params[1]
optimum_alpha = optimum_params[2]

# Plot the optimal solution
if acceptable_params == None:
    print(f"The earliest tumor volume peak occurs after t = {min_peak_time}")
else:

    y0 = y0_base.copy()
    y0[3] = optimum_V0  # Update the initial viral load
    params_modified = POV.params_Ad6
    #params_modified[5] = optimum_e
    params_modified[6] = optimum_c
    params_modified[7] = optimum_alpha
    time = np.linspace(0, 500, num=5000)
    time = np.linspace(0, 500, num=5000)
    earliest_y = odeint(FOV.sim, y0, FOV.t, args=(params_modified,))

    f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)

    LineT, = ax1.plot(FOV.t, earliest_y[:, 0], color="blue", label="LineT")
    LineV, = ax2.plot(FOV.t, earliest_y[:, 3], color="red", label="LineV")

    # plt.plot(t, earliest_y[:, 0], color="blue", label=f"Tumor Volume (V0={earliest_V0:.4f})")
    ax1.set_ylabel('Tumor Volume')
    ax2.set_ylabel('Virus Volume')
    ax2.set_xlabel('Time')
    plt.legend(handles=[LineT, LineV])
    plt.show()
    #print(f"The combination of Optimal V0: {optimum_V0:.4f} and Optimal e: {optimum_e} gives the earliest peak time at {min_peak_time:.2f} days"
       #   f" with a volume of {optimum_peak:.4f} and the tumor remains suppressed for the given time frame.")
    print(
        f"The combination of Optimal V: {optimum_V0:.4f}, Optimal c: {optimum_c}, and Optimal alpha: {optimum_alpha} gives the earliest peak time at {min_peak_time:.2f} days"
        f" with a volume of {optimum_peak:.4f} and the tumor remains suppressed for the given time frame.")

    print(f"The acceptable values of V0 for which the tumor volume peak occurs before t = {FOV.ideal_t_param[-1]}"
          f" and remains suppressed afterwards are {acceptable_params_range}")