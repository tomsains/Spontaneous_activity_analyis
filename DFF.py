import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from oasis.functions import gen_data, gen_sinusoidal_data, deconvolve, estimate_parameters
from oasis.plotting import simpleaxis
from oasis.oasis_methods import oasisAR1, oasisAR2
from scipy import stats
import psutil
import os
import glob



class traces:
    def __init__(self, folder, data_set_name):
        self.folder = folder
        self.data_set_name = data_set_name
        self.traces = np.loadtxt(self.folder + self.data_set_name)
        print("calulating deconvolved")
        self.Apply_DFF()
        self.apply_oasis()


    def apply_oasis(self, pen = 1):
        shape_DFF =  self.DFF.shape
        self.c = np.zeros((shape_DFF))
        self.s = np.zeros((shape_DFF))
        self.c_AR1 = np.zeros((shape_DFF))
        self.s_AR1 = np.zeros((shape_DFF))

        for i in range(shape_DFF[0]):
           self.c_AR1[i, :], self.s_AR1 [i, :] = self.AR1_model_deconvole(self.DFF[i,])
           self.c[i, :], self.s[i, :] = self.deconvolve_trace(self.DFF [i,], penalty=1)

    def deconvolve_trace(self, trace, penalty):
        c, s, b, g, lam = deconvolve(trace, penalty=penalty)
        return (c, s)

    def AR1_model_deconvole(self, trace, smin=0.7):
        c, s = oasisAR1(trace, g=np.exp(-1 / (1.6 * 4.85)), s_min=smin)
        return (c, s)

    def plot_oasis_output(self, cell = range(0,10)):
        fig, axs = plt.subplots(len(cell), 1, figsize=(15, 6), facecolor='w', edgecolor='k')
        fig.subplots_adjust(hspace=.0, wspace=.001)

        axs = axs.ravel()

        for i, c in enumerate(cell):
            print(str(i))
            print(self.b[c])
            axs[i].plot(self.b [c] + self.c [c,:], lw=2, label='denoised')
            axs[i].plot(self.traces [c,:])
            axs[i].plot(self.s [c] - 0.2)
        plt.show()

    def raster_plot(self):
        spikes = self.s > 0.08
        plt.subplot(211)
        plt.imshow(self.c > 0.5, aspect = 2)
        plt.subplot(212)
        plt.plot(np.apply_along_axis(arr= spikes, axis = 0, func1d=sum))



    def base_line(self, t, win = 6000):
        S = pd.Series(t)
        baseline = S.rolling(window = win, center = True, min_periods = 2).quantile(.3)
        return(baseline)




    def DFF_detrend_smooth(self, trace, window = 6000):
        # smooth = butter_lowpass_filter(trace+1, cutoff_freq = 1, nyq_freq = 10/(4.85/2))
        smooth = trace + 1
        baseline = self.base_line(smooth, win=window)
        df = smooth - baseline
        return(df / baseline)

    def Apply_DFF(self):
        self.DFF = np.zeros((self.traces.shape))
        for i in range(self.traces.shape[0]):
            self.DFF[i, :] = np.array(self.DFF_detrend_smooth(self.traces [i,:] ,window=6000))


    def plot_cell(self, number):
        plt.plot(self.DFF [number,:])
        #plt.plot(self.baselines [number,:])
        #plt.plot(self.baselines [number,:] + self.roll_noise [number,:])
        #plt.plot(self.smoothed [number,:])
        #plt.plot(self.noise [number,:])
        plt.show()

    def plot_raster(self):
        plt.imshow(self.DFF)
        plt.show()

    def smoothed_trace(self, t):
        S = pd.Series(t)
        smooth = S.rolling(window=100, win_type='gaussian', center=True).mean(std =5)
        return(smooth)

def get_smooth_DFF(main_folder):
    folders = os.listdir(main_folder)[2:-1]
    for f in folders:
        data_sets = [os.path.basename(x) for x in glob.glob(main_folder + f + "/*traces*")]
        for d in data_sets:
            print("processing ...." + d)
            traces(folder=main_folder + f + "/", data_set_name=d)

def apply_oasis(main_folder, folder = "7_dpf/", result ="/media/thomas_sainsbury/Samsung_T5/SeG/results/Baysian_network_inference/R_scirpts_for_assembly_paper/OASIS_output/WT_NR_7_dpf/", start_from =0 ):
    f = folder
    data_sets = [os.path.basename(x) for x in glob.glob(main_folder + f + "/*traces*")][start_from:]
    for d in data_sets:
            print("processing .... " + d [:-21])
            t = traces(folder=main_folder + f + "/", data_set_name=d)
            np.savetxt(fname=result + d [:-21] + "_oasis_c.txt", X = t.c)
            np.savetxt(fname=result + d[:-21] + "_oasis_s.txt", X = t.s)
            np.savetxt(fname=result + d[:-21] + "_oasisAR1_s.txt", X = t.s_AR1)
            np.savetxt(fname=result + d[:-21] + "_oasisAR1_c.txt", X=t.c_AR1)

