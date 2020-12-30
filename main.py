# This is a sample Python script.

#from remove_intro import *
from DFF import *
from correlations import *
from SNP import *
from synchronous_frames_and_cell_numbers import *
from Utility_functions import load_data_by_deconvolution_method
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import subprocess

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main_folder = "/media/thomas_sainsbury/Samsung_T5/SeG/results/Baysian_network_inference/R_scirpts_for_assembly_paper/BCL_output/"
    decon = False
    AR1_corrs = False
    BCL_corrs = False
    remove_beg = False
    SNPs = False
    sync_frames = False
    test_subampled_corrs = True


    if decon == True:
        apply_oasis(main_folder = main_folder)
        apply_oasis(main_folder = main_folder, folder = "Grav_7_dpf/", result ="/media/thomas_sainsbury/Samsung_T5/SeG/results/Baysian_network_inference/R_scirpts_for_assembly_paper/OASIS_output/WT_GR_7_dpf/")
        apply_oasis(main_folder=main_folder, folder="Grav_5_dpf/",result="/media/thomas_sainsbury/Samsung_T5/SeG/results/Baysian_network_inference/R_scirpts_for_assembly_paper/OASIS_output/WT_GR_5_dpf/")
        apply_oasis(main_folder=main_folder, folder="Grav_3_dpf/",
                    result="/media/thomas_sainsbury/Samsung_T5/SeG/results/Baysian_network_inference/R_scirpts_for_assembly_paper/OASIS_output/WT_GR_3_dpf/")
        apply_oasis(main_folder=main_folder, folder="WT_5_dpf/",
                result="/media/thomas_sainsbury/Samsung_T5/SeG/results/Baysian_network_inference/R_scirpts_for_assembly_paper/OASIS_output/WT_NR_5_dpf/")
        apply_oasis(main_folder=main_folder, folder="WT_3_dpf/",
                   result="/media/thomas_sainsbury/Samsung_T5/SeG/results/Baysian_network_inference/R_scirpts_for_assembly_paper/OASIS_output/WT_NR_3_dpf/")

    if AR1_corrs == True:
        print("calculating_corrs")
        main_folder = "/media/thomas_sainsbury/Samsung_T5/Seg/results/Baysian_network_inference/R_scirpts_for_assembly_paper/OASIS_output/"
        apply_correlations(main_folder = main_folder, folder = "WT_GR_3_dpf")
        apply_correlations(main_folder=main_folder, folder="WT_GR_5_dpf")
        apply_correlations(main_folder=main_folder, folder="WT_GR_7_dpf")
        apply_correlations(main_folder=main_folder, folder="WT_NR_3_dpf")
        apply_correlations(main_folder=main_folder, folder="WT_NR_5_dpf")
        apply_correlations(main_folder=main_folder, folder="WT_GR_7_dpf")

    if BCL_corrs == True:
        main_folder = "/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_output/"
        apply_correlations_all_folder(main_folder=main_folder, deconvolution_method="BCL", iter=140, remove_intro=True, start_from=5)
        apply_correlations_all_folder(main_folder=main_folder, deconvolution_method="BCL", iter=140, remove_intro=False)
        #apply_correlations(main_folder=main_folder, folder="WT_GR_5_dpf/", deconvolution_method="BCL", iter=2, remove_intro= True)
        #apply_correlations(main_folder=main_folder, folder="WT_GR_5_dpf/", deconvolution_method="BCL", iter=2,
            #               remove_intro=False)
        #apply_correlations(main_folder=main_folder, folder="WT_NR_3_dpf/", deconvolution_method="BCL", iter=100, remove_intro=True, start_from=0)
        #apply_correlations(main_folder=main_folder, folder="WT_NR_5_dpf/", deconvolution_method="BCL", iter=140, remove_intro=True, start_from=1)
        #apply_correlations(main_folder=main_folder, folder="WT_NR_7_dpf/", deconvolution_method="BCL", iter=100,
                       #    remove_intro=True, start_from=0)
       # apply_correlations_all_folder(main_folder=main_folder, deconvolution_method="BCL", iter=200, remove_intro=False)
        #apply_correlations_all_folder(main_folder=main_folder, deconvolution_method="AR1", iter=200, remove_intro=True)
        #apply_correlations_all_folder(main_folder=main_folder, deconvolution_method="AR1", iter=200, remove_intro=False)
        #subprocess.call("Rscript  Plot_correlations.R", shell=True)




    if SNPs == True:
        combine_SNPs(main_folder="/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_output/",
                     results_folder="/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_results/SNP/",
                     remove_intro = True, deconvolution_method = "BCL")
        combine_SNPs(main_folder="/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_output/",
                     results_folder="/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_results/SNP/",
                     remove_intro=False, deconvolution_method="BCL")


    if sync_frames == True:
        apply_sync_frames(main_folder="/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_output/",
                          results_folder="/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_results/sync_frames_cell_num/",
                          deconvolution_method="BCL",
                          remove_intro=True)
        apply_sync_frames(main_folder="/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_output/",
                          results_folder="/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_results/sync_frames_cell_num/",
                          deconvolution_method="BCL",
                          remove_intro=False)


    if test_subampled_corrs == True:
        subsample_iter_corrs(folder="/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_output/"+ "WT_GR_7_dpf/", data_set_name="180530_WT_grav_h2b_gc6s_7dpf_f1_sa__00001_scaled_aligned", deconvolution_method="BCL")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
