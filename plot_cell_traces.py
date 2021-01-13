import glob
import os
from Utility_functions import *
import matplotlib.pyplot as plt


class data_set_plot_cell_traces:
    def __init__(self, main_folder, condition_folder, data_set_prefix, deconvolution_method="BCL",
                 remove_intro=True):
        self.main_folder = main_folder
        self.data_set
        self.deconvolution_method
        self.traces = np.load("")
        self.Spikes, self.Calcium = load_data_by_deconvolution_method(full_folder_path=main_folder + condition_folder,
                                                                      data_set_name=data_set_prefix,
                                                                      method=deconvolution_method)

    def plot_cells(self, cell_number=10, start = 500):
        fig, axs = plt.subplots(cell_number, figsize=(20, 15), facecolor='w', edgecolor='w')
        fig.subplots_adjust(hspace=.5, wspace=.001)
        axs = axs.ravel()
        for a, i in enumerate(range(start, start + cell_number)):
            axs[a].plot(self.Calcium[i, :])
            axs[a].plot(self.Spikes[i, :] - 1)
            axs[a].axis("off")

        plt.savefig(self.main_folder + "/cell_traces_" +  )