import numpy as np
import pandas as pd
import

class remove_begining:
    def __init__(self, main_folder, condition_folder, data_set_prefix, deconvolution_method="BCL",
                 remove_intro=" True"):
        if deconvolution_method == "BCL":
            print("loading_cal")
            self.Calcium = np.loadtxt(main_folder + condition_folder + "/" + data_set_prefix + "_all_cells_cal.dat",
                                      dtype="double")
            print("loading_spikes")
            self.Spikes = np.loadtxt(main_folder + condition_folder + "/" + data_set_prefix + "_all_cells_spikes.dat")
            self.Centers = np.loadtxt(main_folder + condition_folder + "/" + data_set_prefix + "_all_cells_spikes.dat")

        self.print_sync_frames_and_cell_number()
        self.remove_intro()
        self.filter_inactive_cells()
        self.print_sync_frames_and_cell_number()
        np.savetxt(fname=main_folder + condition_folder + "/" + data_set_prefix + "_intro_rem_spikes.dat", X= self.Spikes.astype(np.int))
        np.savetxt(fname=main_folder + condition_folder + "/" +  data_set_prefix + "_intro_rem_Cal.dat", X=self.Calcium.astype(np.single))
        np.savetxt(fname=main_folder + condition_folder + "/" + data_set_prefix + "_intro_rem_centers.dat", X=self.Centers.astype(np.int))

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

def remove_intro(main_folder):
    folders = os.listdir(main_folder)[2:-1]

    for f in folders:
        print(f)
        data_set_list = [os.path.basename(x) for x in glob.glob(main_folder + f + "/" + "*_all_cells_spikes.dat")]
        data_sets_prefix = [x.replace("_all_cells_spikes.dat", '') for x in data_set_list]
        for d in data_sets_prefix:
            print("processing ...." + d)
            remove_begining(main_folder=main_folder, condition_folder=f, data_set_prefix=d)