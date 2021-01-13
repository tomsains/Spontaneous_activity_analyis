import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from Utility_functions import load_data_by_deconvolution_method, list_data_sets_in_folder



class sync_frames_cell_num:
    def __init__(self, main_folder, folder, data_set_name,  deconvolution_method, remove_intro, index):
        path = main_folder + folder
        self.spikes, self.calcium = load_data_by_deconvolution_method(full_folder_path=path,data_set_name=data_set_name, method=deconvolution_method)

        if remove_intro == True:
            self.remove_intro()
            self.filter_inactive_cells()


        self.data_set_name = data_set_name
        self.folder = folder
        self.n_active_cells = self.spikes.shape [0]
        self.Rearing_condition = self.folder[3:5]
        self.age = self.folder[6:7]
        self.sync_frames()
        self.build_data_frame(index=index)


    def remove_intro(self, intro=np.int(4.85 * (60 * 5))):
        self.spikes = self.spikes[:, intro:]

    def filter_inactive_cells(self):
        filter_vec = (np.sum(self.spikes, axis=1) > 1)
        self.spikes = self.spikes[filter_vec, :]

    def circular_permutation(self):
        shuff = np.zeros(self.spikes.shape)
        print(shuff.shape)
        for i in range(shuff.shape[0]):
            rand = np.random.uniform(low=0, high=self.spikes.shape[1], size=1)
            shuff[i, :] = np.roll(self.spikes[i, :], int(rand))
        return(shuff)

    def sync_frames(self):
        self.shuff_colSums = np.sum(self.circular_permutation(), axis = 0)
        print(np.quantile(self.shuff_colSums, .95))
        self.real_colSums = np.sum(self.spikes, axis = 0)

    def build_data_frame(self, index):
        self.df = pd.DataFrame({"data_set_name": self.data_set_name, "Rearing_condition": self.Rearing_condition, "age": self.age,
                             "number_of_neurons": self.n_active_cells, "sync_frame_thresh_95": np.quantile(self.shuff_colSums, .95), "sync_frame_thresh_99": np.quantile(self.shuff_colSums,.99),
                             "syc_frames_95": np.sum(self.real_colSums > np.quantile(self.shuff_colSums,.95)),
                             "syc_frames_99": np.sum(self.real_colSums > np.quantile(self.shuff_colSums, .99)),
                             "syc_frames_great_15": np.sum(self.real_colSums > 15)}, index = [index])
        print(self.df ["syc_frames_99"])



def apply_sync_frames(main_folder, results_folder, deconvolution_method, remove_intro):
    print(results_folder + "sync_frame_and_active_cells_" + deconvolution_method + "remove_intro-"+ str(remove_intro) + ".dat")
    folders = [os.path.basename(x) for x in glob.glob(main_folder + "*WT_*")]
    all_conditions = [0]*len(folders)
    for i, f in enumerate(folders):
        print(f)
        data_sets = list_data_sets_in_folder(main_folder=main_folder, condition_folder=f)



        all_fish =  [0]*len(data_sets)
        for j, d in enumerate(data_sets):

            print("processing ...." + d)
            print(d)
            sf_object = sync_frames_cell_num(main_folder=main_folder, folder=f + "/", data_set_name=d,
                                deconvolution_method=deconvolution_method, remove_intro=remove_intro, index=j)

            all_fish [j] = sf_object.df
            #print(all_fish [j])
        all_conditions [i] = pd.concat(all_fish)
    combined_df = pd.concat(all_conditions, ignore_index=True)
    print(combined_df)

    combined_df.to_pickle(results_folder + "sync_frame_and_active_cells_" + deconvolution_method + "_remove_intro-"+ str(remove_intro) + ".pkl")






