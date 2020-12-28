import numpy as np
import pandas as pd

def load_data_by_deconvolution_method(full_folder_path, data_set_name, method):
    if method == "BCL":
        spikes = np.loadtxt(full_folder_path + data_set_name + "_all_cells_spikes.dat")
        calcium = np.loadtxt(full_folder_path + data_set_name + "_all_cells_cal.dat")

    elif method == "AR1":
        spikes = np.loadtxt(full_folder_path + data_set_name + "_oasisAR1_s.txt")
        calcium = np.loadtxt(full_folder_path + data_set_name + "_oasisAR1_c.txt")

    if method == "estimated":
        spikes = np.loadtxt(full_folder_path + data_set_name + "_oasis_s.txt")
        calcium = np.loadtxt(full_folder_path + data_set_name + "_oasis_s.txt")

    return(spikes, calcium)