import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.formula.api import ols
import pingouin as pg
import os
import glob
from Utility_functions import load_data_by_deconvolution_method, list_data_sets_in_folder


class data_frames_for_plots:
    def __init__(self, SNP_folder, sync_frame_folder, deconvolution_method):

        self.sync_frames_intro = pd.read_pickle(sync_frame_folder + "sync_frame_and_active_cells_" + deconvolution_method + "_remove_intro-"+ str(False) + ".pkl")
        self.sync_frames_no_intro = pd.read_pickle(sync_frame_folder + "sync_frame_and_active_cells_" + deconvolution_method + "_remove_intro-" + str(True) + ".pkl")

        self.SNP_intro = pd.read_pickle(SNP_folder + "SNP_table_" + deconvolution_method + "_remove_intro-" + str(False) + ".pkl")
        self.SNP_intro_mean = self.SNP_intro.groupby(["data_set_name", "Rearing_condition", "age"], as_index = False).mean()
        self.SNP_intro_med = self.SNP_intro.groupby(["data_set_name", "Rearing_condition", "age"], as_index=False).median()

        self.SNP_no_intro = pd.read_pickle(SNP_folder + "SNP_table_" + deconvolution_method + "_remove_intro-" + str(True) + ".pkl")
        self.SNP_no_intro_mean = self.SNP_no_intro.groupby(["data_set_name", "Rearing_condition", "age"], as_index=False).mean()
        self.SNP_no_intro_med = self.SNP_no_intro.groupby(["data_set_name", "Rearing_condition", "age"], as_index=False).median()


        self.SNP_folder = SNP_folder
        self.sync_frame_folder = sync_frame_folder
        self.deconvolution_method = deconvolution_method
        self.save_sync_frame_plots()
        self.save_SNP()





    def plot_comparison_sync_frames(self, y, no_intro, intro):
        # sns.set_palette(sns.color_palette(colors))
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))

        colors = {'GR': "#028e96", 'NR': "#ea0037"}
        # sns.boxplot(data = intro, x = "age", y = y, hue = "Rearing_condition", ax =ax[0])
        # sns.boxplot(data = no_intro, x = "age", y = y, hue = "Rearing_condition", ax =ax[1])
        sns.pointplot(x="age", y=y, hue="Rearing_condition", data=intro, ax=ax[0], capsize=0.1, dodge=0.05,
                      palette=colors)
        sns.swarmplot(x="age", y=y, hue="Rearing_condition", data=intro, ax=ax[0], alpha=0.5, palette=colors)
        handles, labels = ax[0].get_legend_handles_labels()
        ax[0].legend(handles[:2], labels[:2])
        ax[0].title.set_text('intro')

        sns.swarmplot(x="age", y=y, hue="Rearing_condition", data=no_intro, ax=ax[1], alpha=0.5, palette=colors)
        sns.pointplot(x="age", y=y, hue="Rearing_condition", data=no_intro, ax=ax[1], capsize=0.1, dodge=0.05,
                      palette=colors, legend=None)
        handles, labels = ax[1].get_legend_handles_labels()
        ax[1].legend(handles[:0], labels[:0])
        ax[1].title.set_text('intro')
        ax[1].title.set_text('no intro')

    def save_sync_frame_plots(self):
        for c in self.sync_frames_intro.columns[3:]:
            print(c)
            self.plot_comparison_sync_frames(y = c, intro=self.sync_frames_intro, no_intro=self.sync_frames_no_intro)
            plt.savefig(self.sync_frame_folder + "plots/sync_fr_and_neuron_numbers_" + self.deconvolution_method + "_" + c + ".pdf")
            plt.close()

            # perform two-way ANOVA
            self.compute_and_save_stats(var = c, data= self.sync_frames_intro, folder = self.sync_frame_folder, intro_removed=False, mean = "mean")
            self.compute_and_save_stats(var=c, data=self.sync_frames_no_intro, folder=self.sync_frame_folder, intro_removed=True, mean="mean")





    def save_SNP(self):
        columns = self.SNP_intro_mean.columns [3:]
    # save_mean_plots
        for c in columns:
            self.plot_comparison_sync_frames(y = c, intro=self.SNP_intro_mean, no_intro=self.SNP_no_intro_mean)
            plt.savefig(self.SNP_folder + "plots/SNP_intro_mean_" + self.deconvolution_method + "_" + c + ".pdf")
            plt.close()

            # save intro stats
            self.compute_and_save_stats(var=c, data=self.SNP_intro_mean, folder=self.SNP_folder, intro_removed=False, mean="mean")
            self.compute_and_save_stats(var=c, data=self.SNP_no_intro_mean, folder=self.SNP_folder, intro_removed=True, mean="mean")


    # save_median_plots
        for c in columns:
            self.plot_comparison_sync_frames(y = c, intro=self.SNP_intro_med, no_intro=self.SNP_no_intro_med)
            plt.savefig(self.SNP_folder + "plots/SNP_intro_med_" + self.deconvolution_method + "_" + c + ".pdf")
            self.compute_and_save_stats(var=c, data=self.SNP_intro_med, folder=self.SNP_folder, intro_removed=False, mean="med")
            self.compute_and_save_stats(var=c, data=self.SNP_no_intro_med, folder=self.SNP_folder, intro_removed=True, mean="med")



    def compute_and_save_stats(self, var, data, folder, intro_removed, mean = "mean"):
        model = ols('{} ~ C(Rearing_condition) + C(age) + C(Rearing_condition):C(age)'.format(var), data=data).fit()
        sm.stats.anova_lm(model, typ=2).to_csv(folder + "stats/intro_removed_" + str(intro_removed) + "_" + self.deconvolution_method + "_" + var + "_" + str(mean) + "_ANOVA.csv")
        t_tests = data.pairwise_ttests(var, between="Rearing_condition", within='age', subject='data_set_name', parametric=True, padjust="bonf", interaction=True)
        t_tests.to_csv(folder + "stats/intro_removed_" + str(intro_removed) + self.deconvolution_method + "_" + var + "_" + str(mean) +  "t_test.csv")




