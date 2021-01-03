import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import os
from scipy import ndimage
from pyrle import Rle
import seaborn as sns
from Utility_functions import load_data_by_deconvolution_method, list_data_sets_in_folder


class data_set_SNP:
    def __init__(self, main_folder, condition_folder, data_set_prefix, deconvolution_method="BCL",
                 remove_intro=" True"):
        self.Spikes, self.Calcium = load_data_by_deconvolution_method(full_folder_path=main_folder+condition_folder,data_set_name=data_set_prefix, method=deconvolution_method)

        self.Centers = np.loadtxt(main_folder + condition_folder + data_set_prefix + "_all_cells_centers_cut_ordered.dat")
        self.print_sync_frames_and_cell_number()
        if remove_intro == True:
            self.remove_intro()
        self.filter_inactive_cells()
        self.print_sync_frames_and_cell_number()

        self.data_set = data_set_prefix
        self.condition_folder = condition_folder
        self.detect_bursts()
        self.sum_active_frames_Hz()
        self.numbered_burst = np.apply_along_axis(arr=self.bursts, axis=1, func1d=self.number_each_burst)
        self.shave_edges_matrix()
        self.burst_stats_matrix()
        self.burst_stats_df = pd.concat([self.burst_stats_df, self.spikes_stats], 1)

    def remove_intro(self, intro=np.int(4.85 * (60 * 5))):
        self.Spikes = self.Spikes[:, intro:]
        self.Calcium = self.Calcium[:, intro:]
        self.Centers = self.Centers[:, intro:]

    def filter_inactive_cells(self):
        filter_vec = (np.sum(self.Spikes, axis=1) > 1)
        self.Spikes = self.Spikes[filter_vec, :]
        self.Calcium = self.Calcium[filter_vec, :]
        self.Centers = self.Centers[filter_vec, :]

    def print_sync_frames_and_cell_number(self):
        # sync_frame = np.sum(self.Spikes, axis = 0)
        print("cell_number:", str(self.Spikes.shape[0]))
        print("Cal_cell_numer:", str(self.Calcium.shape[0]))
        print("sync_frames_above_15_cells:", str(np.sum(np.sum(self.Spikes, axis=0) > 15)))

    def sum_active_frames_Hz(self):
        sum_spikes_Hz = np.sum(self.Spikes, axis=1)
        nonzero_cal_Hz = np.sum(self.Calcium > 0.05, axis=1)
        mean_ISI = np.apply_along_axis(arr=self.Spikes, axis=1, func1d=self.ISI)
        median_ISI = np.apply_along_axis(arr=self.Spikes, axis=1, func1d=self.ISI, stat="median")
        self.spikes_stats = pd.DataFrame(
            {"Spikes_freq": sum_spikes_Hz, "cal_non_zero": nonzero_cal_Hz, "mean_ISI": mean_ISI,
             "median_ISI": median_ISI})

    def plot_cells(self, cell_number, start):
        fig, axs = plt.subplots(cell_number, figsize=(20, 15), facecolor='w', edgecolor='w')
        fig.subplots_adjust(hspace=.5, wspace=.001)
        axs = axs.ravel()
        for a, i in enumerate(range(start, start + cell_number)):
            axs[a].plot(self.Calcium[i, :])
            axs[a].plot(self.Spikes[i, :] - 1)
            axs[a].plot(self.numbered_burst[i,] - 2)
            axs[a].axis("off")

    def ISI(self, trace, stat="mean"):
        a = np.cumsum(trace)
        ISI = np.zeros(int(np.max(a)))
        for i in range(int(np.max(a))):
            ISI[i] = np.sum(a == i)
        if stat == "mean":
            return (np.mean(ISI))
        if stat == "median":
            return (np.median(ISI))

    def detect_bursts(self):
        self.smoothed_spikes = ndimage.gaussian_filter1d(self.Spikes, axis=1, sigma=20)
        self.bursts = self.smoothed_spikes > 0.01

    def number_each_burst(self, trace):
        num_trace = np.zeros(len(trace))
        counter = 1
        for i in range(len(trace) - 1):
            if trace[i] == 1:
                num_trace[i] = counter
            if (trace[i] == 1) & (trace[i + 1] == 0):
                counter = counter + 1
        return (num_trace.astype(np.int))

    def shave_edges(self, numbered_trace, spike_trace):
        new_trace = np.zeros(len(numbered_trace))
        for i in range(1, int(np.max(numbered_trace)) + 1):
            # print(i)
            event = (numbered_trace == i)
            if (sum(spike_trace[event].astype(np.int)) == 0) & (i == int(np.max(numbered_trace)) + 1):
                new_trace[len(spike_trace)] = i

            if (sum(spike_trace[event].astype(np.int)) < 2) & (sum(spike_trace[event].astype(np.int)) > 0):
                place_holder = np.where(numbered_trace == i)[0][0]
                spike_loc = np.where(spike_trace[event])[0][0]

                new_trace[spike_loc + place_holder] = i
                # print(str(spike_loc + place_holder) + "end" )

            if (sum(spike_trace[event].astype(np.int)) >= 2):
                place_holder = np.where(numbered_trace == i)[0][0]
                event_start = np.where(spike_trace[event] == 1)[0][0]
                event_end = np.where(spike_trace[event] == 1)[0][-1]
                new_trace[(place_holder + event_start):(place_holder + event_end)] = i

        return (new_trace)

    def shave_edges_matrix(self):
        for i in range(self.numbered_burst.shape[0]):
            # print(i)
            self.numbered_burst[i, :] = self.shave_edges(numbered_trace=self.numbered_burst[i,],
                                                         spike_trace=self.Spikes[i,])

    def burst_stats(self, numbered_event_trace, spike_trace, cal_trace, index):
        burst_number = np.max(numbered_event_trace.astype(np.int))
        burst_duration = np.zeros(burst_number)
        burst_cal_amp = np.zeros(burst_number)
        burst_sum_spikes = np.zeros(burst_number)
        burst_IBI = self.IBI(numbered_event_trace)
        for i in range(1, burst_number):
            if np.sum(spike_trace[numbered_event_trace == i] > 0) > 0:
                # print(sum(numbered_event_trace == i))
                burst_duration[i] = np.sum(numbered_event_trace == i)
                burst_sum_spikes[i] = np.sum(spike_trace[numbered_event_trace == i])
                burst_cal_amp[i] = np.max(cal_trace[numbered_event_trace == i])
                # print(burst_cal_amp [i])
        return (pd.DataFrame({"data_set_name": [self.data_set], "Rearing_condition": [self.condition_folder[3:5]],
                              "age": [self.condition_folder[6:7]], "burst_number": [burst_number],
                              "burst_duration_mean": [np.mean(burst_duration)],
                              "burst_cal_amp_mean": [np.mean(burst_cal_amp)],
                              "burst_duration_median": np.median(burst_duration),
                              "burst_cal_amp_median": [np.median(burst_cal_amp)],
                              "burst_IBI_mean": [np.mean(burst_IBI)], "burst_IBI_median": [np.median(burst_IBI)]}))

    def burst_stats_matrix(self):
        df_list = [0] * self.Spikes.shape[0]
        for i in range(self.Spikes.shape[0]):
            df_list[i] = self.burst_stats(numbered_event_trace=self.numbered_burst[i,], spike_trace=self.Spikes[i,],
                                          cal_trace=self.Calcium[i,], index=i)
        self.burst_stats_df = pd.concat(df_list, axis=0, ignore_index=True)

    def IBI(self, numbered_trace):
        counts = Rle(numbered_trace)
        return (counts.runs[counts.values < 0.1])


def combine_SNPs(main_folder, results_folder, remove_intro = True, deconvolution_method = "BCL"):
    folders = [os.path.basename(x) for x in glob.glob(main_folder + "*WT_*")]
    print(folders)
    condition_df = [0] * len(folders)
    for j, f in enumerate(folders):
        print(f)
        data_set_list = [os.path.basename(x) for x in glob.glob(main_folder + f + "/" + "*_all_cells_centers_cut_ordered.dat")]
        data_set_prefix = [x.replace("_all_cells_centers_cut_ordered.dat", '') for x in data_set_list]
        print(data_set_prefix)
        data_set_df = [0] * len(data_set_prefix)
        for i, d in enumerate(data_set_prefix):
            print(d)
            d_S = data_set_SNP(main_folder=main_folder, condition_folder=f + "/", data_set_prefix=d, deconvolution_method=deconvolution_method, remove_intro=remove_intro)
            data_set_df[i] = d_S.burst_stats_df
        condition_df[j] = pd.concat(data_set_df)
    total_df = pd.concat(condition_df, ignore_index=True)

    print("saving_at ..... " + results_folder+ "SNP_table_" +  deconvolution_method + "_no_intro.dat")
    total_df.to_pickle(results_folder+ "SNP_table_" +  deconvolution_method + "_remove_intro-" + str(remove_intro) + ".pkl")


