import numpy as np

def find_mean_incidence_angle(peaks, gap=0.4) :
    """
    make a list of the mean incidence angles taken from a list of incidence angles peaks

    Parameters:
    - peaks (array) : array of the incidence angles peaks

    Return: 
    - mean_inc_ang (list) : list of the mean incidence angles taken from the incidence angles peaks list
    """
    peaks = np.concatenate((peaks, [9999]))
    mean_inc_ang = []
    temp_val = []
    for i in range(len(peaks)-1):
        init_peak = peaks[i]
        end_peak = peaks[i+1]
        temp_val.append(init_peak)
        if end_peak - init_peak > gap :     # if the gap between 2 incidence angles is higher than 0.4°
            mean_inc_ang.append(float(np.mean(temp_val)))  # then we mean all the values contained before in temp_val
            temp_val = []                   # we reset the temp_val list   
    mean_inc_ang = np.array(mean_inc_ang)
    return mean_inc_ang