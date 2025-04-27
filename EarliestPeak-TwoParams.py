import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import Parameters_OncolyticVirus as POV
import Functions_OncolyticVirus as FOV

def optimum_and_range_of_V0_and_another(y0_base, params, ideal_t_param, t, param1_range, param2_range):
    min_peak_time = ideal_t_param[-1]
    acceptable_param = [None]
    optimum_params = [None]
    optimum_peak = None
    acceptable_param_range = []

    for param1 in param1_range:
        y0 = y0_base.copy()
        y0[3] = param1  # Update the initial viral load

        for param2 in param2_range:
            params_modified = params.copy()
            params_modified[7] = param2 # Update the virus' sensitivity to interferon


            y = odeint(FOV.sim, y0, t, args=(params_modified,))
            tumor_volume = y[:, 0]

            # Find the peak (maximum) point in tumor volume
            peak_y_index = np.argmax(tumor_volume)
            peak_time = t[peak_y_index]
            peak_y = tumor_volume[peak_y_index]

            # Check if this peak time is the earliest
            if peak_time < ideal_t_param[-1]:
                acceptable_param = [param1, param2]
                acceptable_param_range.append(acceptable_param)
                if peak_time < min_peak_time:
                    min_peak_time = peak_time
                    optimum_params = [param1, param2]
                    optimum_peak = peak_y

    return optimum_params, optimum_peak, acceptable_param_range, acceptable_param, min_peak_time


def optimum_and_range_of_two_nonV0_param(y0_base, params, ideal_t_param, t, param1_range, param2_range):
    min_peak_time = ideal_t_param[-1]
    acceptable_params = [None]
    optimum_params = [None, None]
    optimum_peak = None
    acceptable_param_range = []

    for param1 in param1_range:
        params_modified = params.copy()
        params_modified[6] = param1  # Update the virus' sensitivity to interferon

        for param2 in param2_range:
            params_modified[7] = param2 # Update the virus clearance rate

            y = odeint(FOV.sim, y0_base, t, args=(params_modified,))
            tumor_volume = y[:, 0]

            # Find the peak (maximum) point in tumor volume
            peak_y_index = np.argmax(tumor_volume)
            peak_time = t[peak_y_index]
            peak_y = tumor_volume[peak_y_index]

            # Check if this peak time is the earliest
            if peak_time < ideal_t_param[-1]:
                acceptable_params = [param1, param2]
                acceptable_param_range.append(acceptable_params)
                if peak_time < min_peak_time:
                    min_peak_time = peak_time
                    optimum_params = [param1, param2]
                    optimum_peak = peak_y

    return optimum_params, optimum_peak, acceptable_param_range, acceptable_params, min_peak_time


# Original parameters and settings
y0_base = POV.y0_Ad6  # Base initial conditions with variable V0

# Define a range for V0 and e values
V_range = np.linspace(0.001, 0.02, 20)
e_range = np.linspace(1, 60000, 20)
c_range = np.linspace(0.04, 0.3, 20)
alpha_range = np.linspace(4.854*10**-13, 2.421*10**-7, num=20)

# Find the optimal V0 for the earliest tumor volume peak
optimum_params, optimum_peak, acceptable_params_range, acceptable_params, min_peak_time \
    = optimum_and_range_of_two_nonV0_param(y0_base, POV.params_Ad6,
                                               FOV.ideal_t_param, FOV.t, c_range, alpha_range)


optimum_c = optimum_params[0]
optimum_alpha = optimum_params[1]

# Plot the optimal solution
if acceptable_params == None:
    print(f"The earliest tumor volume peak occurs after t = {min_peak_time}")
else:

    y0 = y0_base.copy()
    #y0[3] = optimum_V0  # Update the initial viral load
    params_modified = POV.params_Ad6
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
    print(
        f"The combination of Optimal c: {optimum_c:.4f}, Optimal alpha: {optimum_alpha} gives the earliest peak time at {min_peak_time:.2f} days"
        f" with a volume of {optimum_peak:.4f} mm^3 and the tumor remains suppressed for the given time frame.")

    print(f"The acceptable values of V0 for which the tumor volume peak occurs before t = {FOV.ideal_t_param[-1]}"
          f" and remains suppressed afterwards are {acceptable_params_range}")