import numpy as np
import pandas as pd
import glob
import os


def load_data_by_deconvolution_method(full_folder_path, data_set_name, method):
    if method == "BCL":
        spikes = np.loadtxt(full_folder_path + data_set_name + "_all_cells_spikes.dat")
        calcium = np.loadtxt(full_folder_path + data_set_name + "_all_cells_cal.dat")

    elif method == "AR1":
        spikes = np.loadtxt(full_folder_path + data_set_name + "_oasisAR1_s.txt")
        calcium = np.loadtxt(full_folder_path + data_set_name + "_oasisAR1_c.txt")
        spikes = (spikes > 0.4) * 1

    elif method == "estimated":
        spikes = np.loadtxt(full_folder_path + data_set_name + "_oasis_s.txt")
        calcium = np.loadtxt(full_folder_path + data_set_name + "_oasis_s.txt")
        spikes = (spikes > 0.1) * 1

    return(spikes, calcium)

def list_data_sets_in_folder(main_folder, condition_folder):
    data_set_list = [os.path.basename(x) for x in glob.glob(main_folder + condition_folder + "/" + "*_all_cells_centers_cut_ordered.dat")]
    data_set_prefix = [x.replace("_all_cells_centers_cut_ordered.dat", '') for x in data_set_list]
    return(data_set_prefix)