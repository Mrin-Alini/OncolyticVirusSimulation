import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the differential equations
def sim(variables, t, params):
    T, E, I, V, F = variables
    lam, beta, k, ro, p, e, c, alpha = params

    dTdt = lam * T - beta * T * V
    dEdt = beta * T * V - k * E
    dIdt = k * E - ro * I
    dVdt = I * (p / (1 + e * F)) - c * V
    dFdt = V - alpha * F

    return [dTdt, dEdt, dIdt, dVdt, dFdt]

# Time points with increased resolution
t = np.linspace(0, 365, num=3650)
ideal_t_param = np.linspace(0, 365, num=3650)

ideal_v_param = np.linspace(0, 3000, num=30000)

# Finding Earliest Peak
def find_earliest_peak(y0_base, params, ideal_t_param, param_range):
    min_peak_time = ideal_t_param[-1]
    earliest_param = None
    peak_y = None

    for param in param_range:
        y0 = y0_base.copy()
        y0[3] = param  # Update the initial viral load

        y = odeint(sim, y0, ideal_t_param, args=(params,))
        tumor_volume = y[:, 0]

        # Find the peak (maximum) point in tumor volume
        peak_y_index = np.argmax(tumor_volume)
        peak_time = ideal_t_param[peak_y_index]

        # Check if this peak time is the earliest
        if peak_time < min_peak_time:
            min_peak_time = peak_time
            earliest_param = param
            peak_y = tumor_volume[peak_y_index]

    return earliest_param, min_peak_time, peak_y

# Finding Earliest Peak that is Suppressed
def find_earliest_peak_and_suppression(y0_base, params, ideal_t_param, t, param_range):
    min_peak_time = ideal_t_param[-1]
    earliest_param = None
    peak_y = None

    for param in param_range:
        y0 = y0_base.copy()
        y0[3] = param # Update the initial viral load
        y = odeint(sim, y0, ideal_t_param, args=(params,))
        tumor_volume = y[:, 0]

        # Find the peak (maximum) point in tumor volume
        peak_y_index = np.argmax(tumor_volume)
        peak_time = ideal_t_param[peak_y_index]

        # Check if this peak time is the earliest
        if peak_time < min_peak_time:
            min_peak_time = peak_time
            earliest_param = param
            peak_y = tumor_volume[peak_y_index]

    y0 = y0_base.copy()
    y0[3] = earliest_param
    Y = odeint(sim, y0, t, args=(params,))
    total_tumor_growth = Y[:, 0]

    if max(total_tumor_growth) == peak_y:
        suppression = True
    else:
        suppression = False

    return earliest_param, min_peak_time, peak_y, suppression, total_tumor_growth


# Finding Optimum and Acceptable Parameters
#def find_optimum_and_range_of_param(y0_base, params, ideal_t_param, t, param_range):
    #min_peak_time = ideal_t_param[-1]
    #acceptable_param = None
    #optimum_param = None
    #optimum_peak = None
    #acceptable_param_range = []

    #for param in param_range:
        #y0 = y0_base.copy()
        #y0[3] = param  # Update the initial viral load
        #y = odeint(sim, y0, ideal_t_param, args=(params,))
        #tumor_volume = y[:, 0]

        # Find the peak (maximum) point in tumor volume
        #peak_y_index = np.argmax(tumor_volume)
        #peak_time = ideal_t_param[peak_y_index]
        #peak_y = tumor_volume[peak_y_index]

        # Check if this peak time is the earliest
        #if peak_time < ideal_t_param[-1]:
            #acceptable_param = param
            #Y = odeint(sim, y0, t, args=(params,))
            #total_tumor_growth = Y[:, 0]
            #total_largest_index = np.argmax(total_tumor_growth)
            #total_largest = total_tumor_growth[total_largest_index]
            #print(param, peak_time, ideal_t_param[-1], peak_y, total_largest)
            #if abs((total_largest - peak_y)) < 0.001 :
                #acceptable_param_range.append(acceptable_param)
                #if peak_time < min_peak_time:
                    #min_peak_time = peak_time
                    #optimum_param = param
                    #optimum_peak = peak_y


    #return optimum_param, optimum_peak, acceptable_param_range, acceptable_param, min_peak_time




def find_optimum_and_range_of_param(y0_base, params, ideal_t_param, t, param_range):
    min_peak_time = ideal_t_param[-1]
    acceptable_param = None
    optimum_param = None
    optimum_peak = None
    acceptable_param_range = []

    for param in param_range:
        y0 = y0_base.copy()
        param_modified = params.copy()
        param_modified[7] = param # Update the alpha
        y = odeint(sim, y0, t, args=(param_modified,))
        tumor_volume = y[:, 0]

        # Find the peak (maximum) point in tumor volume
        peak_y_index = np.argmax(tumor_volume)
        peak_time = t[peak_y_index]
        peak_y = tumor_volume[peak_y_index]

        # Check if this peak time is the earliest
        if peak_time < ideal_t_param[-1]:
            acceptable_param = param
            acceptable_param_range.append(acceptable_param)
            if peak_time < min_peak_time:
                min_peak_time = peak_time
                optimum_param = param
                optimum_peak = peak_y


    return optimum_param, optimum_peak, acceptable_param_range, acceptable_param, min_peak_time


def find_optimum_and_range_of_many_param(y0_base, params, ideal_t_param, t, param1_range, param2_range, param3_range):
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
            params_modified[6] = param2 # Update the virus' sensitivity to interferon

            for param3 in param3_range:
                params_modified[7] = param3

                y = odeint(sim, y0, t, args=(params_modified,))
                tumor_volume = y[:, 0]

                # Find the peak (maximum) point in tumor volume
                peak_y_index = np.argmax(tumor_volume)
                peak_time = t[peak_y_index]
                peak_y = tumor_volume[peak_y_index]

                # Check if this peak time is the earliest
                if peak_time < ideal_t_param[-1]:
                    acceptable_param = [param1, param2, param3]
                    acceptable_param_range.append(acceptable_param)
                    if peak_time < min_peak_time:
                        min_peak_time = peak_time
                        optimum_params = [param1, param2, param3]
                        optimum_peak = peak_y

    return optimum_params, optimum_peak, acceptable_param_range, acceptable_param, min_peak_time


def find_optimum_and_range_of_multiple_param(y0_base, params, ideal_t_param, t, param1_range, param2_range, param3_range):
    min_peak_time = ideal_t_param[-1]
    acceptable_params = [None]
    optimum_params = [None, None]
    optimum_peak = None
    acceptable_param_range = []

    for param1 in param1_range:
        params_modified = params.copy()
        params_modified[5] = param1  # Update the virus' sensitivity to interferon

        for param2 in param2_range:
            params_modified[6] = param2 # Update the virus clearance rate

            for param3 in param3_range:
                params_modified[7] = param2  # Update the alpha

                y = odeint(sim, y0_base, t, args=(params_modified,))
                tumor_volume = y[:, 0]

              # Find the peak (maximum) point in tumor volume
                peak_y_index = np.argmax(tumor_volume)
                peak_time = t[peak_y_index]
                peak_y = tumor_volume[peak_y_index]

            # Check if this peak time is the earliest
                if peak_time < ideal_t_param[-1]:
                    acceptable_params = [param1, param2, param3]
                    acceptable_param_range.append(acceptable_params)
                    if peak_time < min_peak_time:
                        min_peak_time = peak_time
                        optimum_params = [param1, param2, param3]
                        optimum_peak = peak_y

    return optimum_params, optimum_peak, acceptable_param_range, acceptable_params, min_peak_time


def find_optimum_and_range_of_all_four(y0_base, params, ideal_t_param, t, param1_range, param2_range, param3_range, param4_range):
    min_peak_time = ideal_t_param[-1]
    acceptable_params = [None]
    optimum_params = [None]
    optimum_peak = None
    acceptable_param_range = []

    for param1 in param1_range:
        y0 = y0_base.copy()
        y0[3] = param1  # Update the initial viral load

        for param2 in param2_range:
            params_modified = params.copy()
            params_modified[5] = param2 # Update the virus' sensitivity to interferon

            for param3 in param3_range:
                params_modified[6] = param3 # Update the virus clearance rate

                for param4 in param4_range:
                    params_modified[7] = param4 # Update decay rate of interferon

            y = odeint(sim, y0, t, args=(params_modified,))
            tumor_volume = y[:, 0]

            # Find the peak (maximum) point in tumor volume
            peak_y_index = np.argmax(tumor_volume)
            peak_time = t[peak_y_index]
            peak_y = tumor_volume[peak_y_index]

            # Check if this peak time is the earliest
            if peak_time < ideal_t_param[-1]:
                acceptable_params = [param1, param2, param3, param4]
                acceptable_param_range.append(acceptable_params)
                if peak_time < min_peak_time:
                    min_peak_time = peak_time
                    optimum_params = [param1, param2, param3, param4]
                    optimum_peak = peak_y

    return optimum_params, optimum_peak, acceptable_param_range, acceptable_params, min_peak_time



def find_optimum_and_range_of_all_four_volume(y0_base, params, ideal_v_param, t, param1_range, param2_range, param3_range, param4_range):
    min_volume = ideal_v_param[-1]
    acceptable_params = [None]
    optimum_params = [None]
    optimum_peak_time = None
    acceptable_param_range = []

    for param1 in param1_range:
        y0 = y0_base.copy()
        y0[3] = param1  # Update the initial viral load

        for param2 in param2_range:
            params_modified = params.copy()
            params_modified[5] = param2 # Update the virus' sensitivity to interferon

            for param3 in param3_range:
                params_modified[6] = param3 # Update the virus clearance rate

                for param4 in param4_range:
                    params_modified[7] = param4 # Update decay rate of interferon

            y = odeint(sim, y0, t, args=(params_modified,))
            tumor_volume = y[:, 0]

            # Find the peak (maximum) point in tumor volume
            peak_y_index = np.argmax(tumor_volume)
            peak_time = t[peak_y_index]
            peak_y = tumor_volume[peak_y_index]

            # Check if this peak time is the earliest
            if peak_y < ideal_v_param[-1]:
                acceptable_params = [param1, param2, param3, param4]
                acceptable_param_range.append(acceptable_params)
                if peak_y < min_volume:
                    min_volume = peak_y
                    optimum_params = [param1, param2, param3, param4]
                    optimum_peak_time = peak_time

    return optimum_params, min_volume, acceptable_param_range, acceptable_params, optimum_peak_time
