import numpy as np
import pandas as pd
import psutil
import os
import glob
from sklearn.metrics.pairwise import euclidean_distances
import matplotlib.pyplot as plt
from Utility_functions import load_data_by_deconvolution_method


class correlations:
    def __init__(self, folder, data_set_name, deconvolution_method="AR1", iter = 200, remove_intro =True, midline_no = 0):
        self.spikes, self.calcium = load_data_by_deconvolution_method(full_folder_path=folder,data_set_name=data_set_name, method=deconvolution_method)

        if remove_intro == True:
            self.remove_intro()

        self.centers = np.loadtxt(folder + data_set_name + "_all_cells_centers_cut_ordered.dat")[:,0:2]

        self.filter_inactive_cells()
        print("centers shape:" + str(self.centers.shape))
        self.midline = pd.read_csv(folder + "midlines.dat", delimiter= " ").iloc [midline_no]
        print(self.midline)
        print("spikes shape:" + str(self.spikes.shape))
        #self.plot_centers_with_midline()
        self.null(iter=iter)

    def filter_inactive_cells(self):
        filter_vec = (np.sum(self.spikes, axis=1) > 1)
        self.spikes = self.spikes[filter_vec, :]
        self.calcium = self.calcium[filter_vec, :]
        self.centers = self.centers[filter_vec, :]

    def circular_permutation(self):
        shuff = self.spikes
        for i in range(shuff.shape[0]):
            rand = np.random.uniform(low=0, high=self.spikes.shape[0], size=1)
            shuff[i, :] = np.roll(shuff[i, :], int(rand))
        return (shuff)

    def plot_centers_with_midline(self):
        plt.scatter(self.centers[:,0], self.centers[:,1])
        abline_values = [self.midline ["slope"] * i + self.midline ["intercept"] for i in self.centers[:,0]]
        plt.plot(self.centers [:,0], abline_values)
        plt.show()

    def which_side_is_cell(self):
        hemisphere = np.zeros(self.centers.shape [0])
        hemisphere [self.centers [:,1] > (self.centers [:,0]*self.midline ["slope"] + self.midline ["intercept"])] = 1
        hemisphere[self.centers[:, 1] < (self.centers[:, 0] * self.midline ["slope"] + self.midline ["intercept"])] = -1
        self.side = hemisphere
        print(sum(self.side == 1))
        print(sum(self.side == -1))
        #plt.scatter(self.centers [self.side == 1,0],self.centers [self.side == 1,1])
        #plt.scatter(self.centers [self.side == -1,0],self.centers [self.side == -1,1])
        #plt.plot(self.centers [:,0], self.centers [:,0]*self.midline ["slope"] + self.midline ["intercept"], c = "red")
        #plt.show()

    def null(self, iter=2):
        self.nullcorrs = np.zeros(shape=(self.spikes.shape[0], self.spikes.shape[0], iter))
        self.realcorrs = np.corrcoef(self.spikes)
        self.calcorrs = np.corrcoef(self.calcium)
        self.cell_distances = euclidean_distances(self.centers)
        self.which_side_is_cell()
        print("real_corrs_shape" + str(self.realcorrs.shape))
        print(self.cell_distances.shape)

        for i in range(iter):
            print(i)
            shuff = self.circular_permutation()
            self.nullcorrs[:, :, i] = np.corrcoef(shuff)
            if i % 10 == 0:
                print(psutil.virtual_memory())
        print("calculating_quantile")
        self.nullcorrs = np.apply_along_axis(func1d=np.quantile, arr=self.nullcorrs, axis=2, q=.99)

    def remove_intro(self, intro=np.int(4.85 * (60 * 5))):
        self.spikes = self.spikes[:, intro:]

'''
def get_corrs(main_folder):
    folders = os.listdir(main_folder)[2:-1]
    for f in folders:
        data_sets = [os.path.basename(x) for x in glob.glob(main_folder + f + "/*traces*")]
        for d in data_sets:
            print("processing ...." + d)
            traces(folder=main_folder + f + "/", data_set_name=d)
'''

def apply_correlations(main_folder, folder="WT_GR_3_dpf",deconvolution_method = "AR1", iter = 120, remove_intro = True, start_from = 0):
        f = folder
        if deconvolution_method == "AR1":
            data_sets = [os.path.basename(x) for x in glob.glob(main_folder + f + "/*_oasis_c*")]
            suffix_len = 12

        if deconvolution_method == "BCL":
            data_sets = [os.path.basename(x) for x in glob.glob(main_folder + f + "/*_spikes*")]
            print(data_sets)

            suffix_len = 22 -1

        print(len(data_sets))
        data_sets = data_sets [start_from:len(data_sets)]
        for i, d in enumerate(data_sets):

            print("processing ...." + d[:-suffix_len])
            print(d)
            corrs = correlations(folder=main_folder + f, data_set_name= d [:-suffix_len], deconvolution_method = deconvolution_method, iter=iter, remove_intro=remove_intro, midline_no=i)
            print("save_real_corrs")
            np.savetxt(fname=main_folder + folder + "/correlations/" + d[:-suffix_len] + deconvolution_method + "_cal_real_correlations_removed_intro-" + str(remove_intro) + ".npy", X=corrs.calcorrs)
            np.savetxt(fname=main_folder + folder + "/correlations/" + d[:-suffix_len] + deconvolution_method + "_spikes_real_correlations_removed_intro-" + str(remove_intro) + ".npy", X=corrs.realcorrs)

            print("save_null_corrs")
            np.savetxt(fname=main_folder + folder + "/correlations/" + d[:-suffix_len] + deconvolution_method + "_spikes_null_correlations_sig_0.99_iter_" + str(iter) + "_removed_intro-" + str(remove_intro) + ".npy", X=corrs.nullcorrs)

            print("save_distances")
            np.savetxt(fname=main_folder + folder + "/correlations/" + d[:-suffix_len] + deconvolution_method + "_cell_distances" + "_removed_intro-" + str(remove_intro) + ".npy", X=corrs.cell_distances)

            print()
            np.savetxt(fname=main_folder + folder + "/correlations/" + d[:-suffix_len] + deconvolution_method + "_cell_side" + "_removed_intro-" + str(remove_intro) + ".npy", X=corrs.side.astype(np.int))
            del corrs


def apply_correlations_all_folder(main_folder, deconvolution_method, iter, remove_intro):
    folders = [os.path.basename(x) for x in glob.glob(main_folder + "*WT_*")]
    for f in folders:
        print(f)
        apply_correlations(main_folder=main_folder, folder=f +"/", deconvolution_method = deconvolution_method, iter = iter, remove_intro=remove_intro)
