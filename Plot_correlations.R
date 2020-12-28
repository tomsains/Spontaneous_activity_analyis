# Title     : TODO
# Objective : TODO
# Created by: thomas_sainsbury
# Created on: 21/12/2020

library(data.table)
#library(sm)
library(ggplot2)
library(viridis)
library(DescTools)
library(fields)
library(matrixStats)
library(plyr)
library(data.table)
library(sm)
library(ggplot2)
library(viridis)
library(DescTools)
library(fields)
library(matrixStats)
library(car)
library(Hmisc)
library(rdist)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(GCalcium)
library(foreach)
library(stringr)
library(doParallel)

setwd("/media/thomas_sainsbury/Samsung_T5/Spontaneous_activity_experiments/BCL_output/")

load_corr <- function(folder = "WT_GR_7_dpf/correlations/", data_set = "180530_WT_grav_h2b_gc6s_7dpf_f1_sa__00001_scaled_aligned", deconvolution_method = "BCL", remove_intro = "True", iter = 140) {
  corr_cal <- as.matrix(fread(paste(folder, data_set, deconvolution_method, "_cal_real_correlations_removed_intro-", remove_intro, ".dat", sep ="")))
  corr_spikes <- as.matrix(fread(paste(folder, data_set, deconvolution_method, "_spikes_real_correlations_removed_intro-", remove_intro, ".dat", sep ="")))
  side <- as.matrix(fread(paste(folder, data_set, deconvolution_method, "_cell_side", "_removed_intro-", remove_intro, ".dat", sep ="")))
  NULL_corrs <- as.matrix(fread(paste(folder, data_set, deconvolution_method, "_spikes_null_correlations_sig_0.99_iter_", iter, "_removed_intro-", remove_intro, ".dat", sep ="")))
  distance <- as.matrix(fread(paste(folder, data_set, deconvolution_method, "_cell_distances", "_removed_intro-" , remove_intro,".dat", sep ="")))

  same_side <- side %*% t(side)
  return(list("spikes_corr" = corr_spikes, "cal_corr" = corr_cal, "null_corrs" = NULL_corrs, "distance_between_cells" = distance*1.3, "same_side" = same_side,  "side" = side))
}

get_lower_tri_mat <- function(mat){
  return(mat[lower.tri(mat)])
}

get_lowers_tri_for_corrs <- function(corrs =corrs) {
  corrs$spikes_corr <- get_lower_tri_mat(corrs$spikes_corr)
  corrs$cal_corr <- get_lower_tri_mat(corrs$cal_corr)
  corrs$null_corrs <- get_lower_tri_mat(corrs$null_corrs)
  corrs$distance_between_cells <- get_lower_tri_mat(corrs$distance_between_cells)
  corrs$same_side <- get_lower_tri_mat( corrs$same_side)
  return(corrs)
}

cut_corrs_by_dist <- function(corrs = corrs){

  same_side_cor_spikes <- by(corrs$spikes_corr [(corrs$spikes > corrs$null_corrs) & (corrs$same_side == 1)], cut(corrs$distance_between_cells[(corrs$spikes > corrs$null_corrs) & (corrs$same_side == 1)]*1.3, breaks =seq(0, 300, 50)), mean, na.rm = TRUE)
  same_side_cor_cal <- by(corrs$cal_corr [(corrs$spikes > corrs$null_corrs) & (corrs$same_side == 1)], cut(corrs$distance_between_cells[(corrs$spikes > corrs$null_corrs) & (corrs$same_side == 1)]*1.3, breaks =seq(0, 300, 50)), mean, na.rm = TRUE)
  return(c(as.vector(same_side_cor_spikes), as.vector(same_side_cor_cal)))
}


get_lower_tri_mat <- function(mat){
  return(mat[lower.tri(mat)])
}


pairwise_correlation_summaries <- function(folder = "WT_GR_7_dpf/correlations/", data_set = "180530_WT_grav_h2b_gc6s_7dpf_f1_sa__00001_scaled_aligned", remove_intro){
  corrs <- load_corr(folder  = folder, data_set = data_set, remove_intro = remove_intro)
  data <- get_lowers_tri_for_corrs(corrs)
  dat <- cbind(mean(data$spikes_corr [(data$spikes > data$null_corrs) & (data$same_side == 1)], na.rm = TRUE),
               mean(data$cal_corr [(data$spikes > data$null_corrs) & (data$same_side == 1)], na.rm = TRUE),
               mean(data$spikes_corr [(data$spikes > data$null_corrs) & (data$same_side == -1)], na.rm = TRUE),
               mean(data$cal_corr [(data$spikes > data$null_corrs) & (data$same_side == -1)], na.rm = TRUE),
                median(data$spikes_corr [(data$spikes > data$null_corrs) & (data$same_side == 1)], na.rm = TRUE),
                median(data$cal_corr [(data$spikes > data$null_corrs) & (data$same_side == 1)], na.rm = TRUE),
                median(data$spikes_corr [(data$spikes > data$null_corrs) & (data$same_side == -1)], na.rm = TRUE),
                median(data$cal_corr [(data$spikes > data$null_corrs) & (data$same_side == -1)], na.rm = TRUE),
               t(cut_corrs_by_dist(data)))
   colnames(dat) <- c("mean_corr_spikes_ipsi", "mean_corr_calcium_ipsi", "mean_corr_spikes_contra", "mean_corr_calcium_contra",
                     "median_corr_spikes_ipsi", "median_corr_calcium_ipsi", "median_corr_spikes_conta", "median_corr_spikes_contra",
                     sprintf("Distance_spike_corr_%03d", seq(0,250,50)),sprintf("Distance_cal_corr_%03d", seq(0,250,50)))

  return(dat)
}

list_all_fish_in_dir <- function(folder, deconvolution_method = "BCL", remove_intro ){
  list_fish = list.files(path = folder, pattern = paste(deconvolution_method, "_cell_side_removed_intro-", remove_intro, ".dat", sep =""))
  return(lapply(list_fish, str_replace, pattern = paste(deconvolution_method, "_cell_side_removed_intro-", remove_intro, ".dat", sep =""), replacement = ""))
}

corrs_all_fish_in_group <- function(folder, genotype, age, rearing_condition, deconvolution_method = "BCL", remove_intro = "True"){
  folder <- paste(folder, "correlations/", sep ="")
  data_sets <- list_all_fish_in_dir(folder=folder, deconvolution_method, remove_intro = remove_intro)
  print(data_sets)
  cor_for_group <- vector("list", length(data_sets))
  for (i in 1:length(data_sets)){
    cor_for_group [[i]] <- pairwise_correlation_summaries(folder = folder, data_set = data_sets[i], remove_intro = remove_intro)
  }
  cor_for_group <- as.data.frame(do.call(rbind, cor_for_group))
  group_stats <- data.frame("data_set_names" = unlist(data_sets), "Genotype" = rep(genotype, length(data_sets)), "Age" = rep(age, length(data_sets)), "Rearing_condition" = rep(rearing_condition, length(data_sets)))
  dat <-cbind(group_stats,cor_for_group)

  return(dat)
}




plot_correlations_for_each_age <- function(age = "5_dpf", save_folder = "../BCL_results/Correlation_results/BCL_intro_removed/") {

  if (age == "3_dpf") {
    combine <- rbind(o_norm_3, o_grav_3)
  }
  if (age == "5_dpf") {
    combine <- rbind(o_norm_5, o_grav_5)
  }
  if (age == "7_dpf") {
    combine <- rbind(o_norm_7, o_grav_7)
  }
  combine$ipsi_contra_ratio <- combine$mean_corr_spikes_ipsi/combine$mean_corr_spikes_contra
  combine$mean_corr_spikes_distance_difference <- (c(combine$Distance_spike_corr_000 + combine$Distance_spike_corr_100 + combine$Distance_spike_corr_050)/3)/(c(combine$Distance_spike_corr_250 + combine$Distance_spike_corr_200)/2)


  combine %>% select(Rearing_condition, sprintf("Distance_spike_corr_%03d", seq(0,250,50))) %>%
    pivot_longer(., cols = c(sprintf("Distance_spike_corr_%03d", seq(0,250,50))), names_to = "Var", values_to = "Val") %>%
    ggplot(aes(x = Var, y = Val, fill = Rearing_condition)) +
    geom_boxplot() + xlab("Distance (um)") + ylab("Correlation Coef.") + scale_x_discrete(labels=c("0-50", "50-100", "100-150", "150-200", "200-250", "250-300"))

  ggsave(paste(save_folder, "Distance_spikes_corrs_", age, ".pdf", sep = ""))

  combine %>% select(Rearing_condition, sprintf("Distance_cal_corr_%03d", seq(0,250,50))) %>%
    pivot_longer(., cols = c(sprintf("Distance_cal_corr_%03d", seq(0,250,50))), names_to = "Var", values_to = "Val") %>%
    ggplot(aes(x = Var, y = Val, fill = Rearing_condition)) +
    geom_boxplot() + xlab("Distance (um)") + ylab("Correlation Coef.") + scale_x_discrete(labels=c("0-50", "50-100", "100-150", "150-200", "200-250", "250-300"))
  ggsave(paste(save_folder, "Distance_cal_corrs_", age, ".pdf", sep = ""))

  combine %>% select(Rearing_condition, mean_corr_spikes_ipsi, mean_corr_spikes_contra) %>%
    pivot_longer(., cols = c(mean_corr_spikes_ipsi, mean_corr_spikes_contra), names_to = "Var", values_to = "Val") %>%
    ggplot(aes(x = Var, y = Val, fill = Rearing_condition)) +
    geom_boxplot()  + xlab("Tectal hemishpere") + ylab("Correlation Coef.") + scale_x_discrete(labels=c("Contra", "Ipsi"))
  ggsave(paste(save_folder,"Mean_spike_corrs_ipsi_vs_contra_", age, ".pdf", sep = ""))

  combine %>% select(Rearing_condition, mean_corr_spikes_distance_difference) %>%
    pivot_longer(., cols = c(mean_corr_spikes_distance_difference), names_to = "Var", values_to = "Val") %>%
    ggplot(aes(x = Var, y = Val, fill = Rearing_condition)) +
    geom_boxplot() + xlab("Diff. Short vs long range Correlations") + ylab("Correlation Coef.")
  ggsave(paste(save_folder, "Diff_Short_vs_long_range_corr_", age, ".pdf", sep = ""))

  combine %>% select(Rearing_condition, ipsi_contra_ratio) %>%
    pivot_longer(., cols = c(ipsi_contra_ratio), names_to = "Var", values_to = "Val") %>%
    ggplot(aes(x = Var, y = Val, fill = Rearing_condition)) +
    geom_boxplot() + xlab("IPSI_vs_contra_ratio") + ylab("Correlation Coef.")
  ggsave(paste(save_folder, "IPSI_vs_contra_Correlation_ratio_", age, ".pdf", sep = ""))
}

plot_correltaions_by_age <- function(results_folder){
  plot_correlations_for_each_age(age = "3_dpf", save_folder = results_folder)
  plot_correlations_for_each_age(age = "5_dpf", save_folder  = results_folder)
  plot_correlations_for_each_age(age = "7_dpf", save_folder = results_folder)
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(as.numeric(x[[col]]), na.rm=TRUE))
    c(sd = sd(as.numeric(x[[col]]), na.rm=TRUE))
    c(med = median(as.numeric(x[[col]]), na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

plot_line_corrs_across_age <-function(data = dev , y="Distance_spike_corr_000", ylab = "Short range Corr coeff",  col = c("#ea0037","#028e96")){
  data %>% select("Age", "Rearing_condition", all_of(y)) %>% ggplot(aes_string(x = "Age", y= y, group = "Rearing_condition", col = "Rearing_condition"))  + ylab(ylab) + ylim(ylim)  + stat_summary(fun=mean, geom="line", lwd = 1)  +
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.2,lwd = 1) +
    stat_summary(fun=mean, geom="point", size = 5) + xlab("Age (dpf)") + theme_classic(base_size = 40) +
    scale_color_manual(values = col) + theme(legend.position="none") + geom_jitter(width = 0.05, alpha =0.5, size = 3)
}

plot_corrs_across_age <- function(save_folder = "../BCL_results/Correlation_results/BCL_intro_removed/"){
  dev <- as.data.frame(rbind(o_norm_3, o_grav_3, o_norm_5, o_norm_5, o_grav_5, o_norm_7, o_grav_7))
  dev$short_range_correlations <- apply(cbind(dev$Distance_spike_corr_000, dev$Distance_spike_corr_100, dev$Distance_spike_corr_150),1,mean)
  dev$long_range_correlations <- apply(cbind(dev$Distance_spike_corr_250, dev$Distance_spike_corr_200),1,mean)
  dev$long_short_ratio_spikes <- dev$short_range_correlations/dev$long_range_correlations

  plot_line_corrs_across_age(data = dev, y = "short_range_correlations")
  ggsave(paste(save_folder, "Short_range_corrs_across_development", ".pdf", sep = ""))

  plot_line_corrs_across_age(data = dev,y = "long_range_correlations", ylab = "Long range corr coeff")
  ggsave(paste(save_folder, "Long_range_corrs_across_development", ".pdf", sep = ""))

  plot_line_corrs_across_age(data = dev,y = "long_short_ratio_spikes", ylab = "short/long range corrs")
  ggsave(paste(save_folder, "Long_short_corrs_ratio_across_development", ".pdf", sep = ""))

  dev$short_range_correlations_cal <- apply(cbind(dev$Distance_cal_corr_000, dev$Distance_cal_corr_050),1,mean)
  dev$long_range_correlations_cal <- apply(cbind(dev$Distance_cal_corr_200, dev$Distance_cal_corr_250),1,mean)
  dev$long_short_ratio_cal <- dev$short_range_correlations_cal/dev$long_range_correlations_cal

  plot_line_corrs_across_age(data = dev, y = "short_range_correlations_cal")
  ggsave(paste(save_folder, "Short_range_corrs_across_development_cal", ".pdf", sep = ""))

  plot_line_corrs_across_age(data = dev,y = "long_range_correlations_cal", ylab = "Long range corr coeff")
  ggsave(paste(save_folder, "Long_range_corrs_across_development_cal", ".pdf", sep = ""))

  plot_line_corrs_across_age(data = dev,y = "long_short_ratio_cal", ylab = "short/long range corrs")
  ggsave(paste(save_folder, "Long_short_corrs_ratio_across_development_cal", ".pdf", sep = ""))

  plot_line_corrs_across_age(data = dev,y = "mean_corr_spikes_ipsi", ylab = "mean corr Coeff. (spikes)")
  ggsave(paste(save_folder, "mean_corr_spikes_ipsi", ".pdf", sep = ""))

  plot_line_corrs_across_age(data = dev,y = "mean_corr_calcium_ipsi", ylab = "mean corr Coeff. (cal)")
  ggsave(paste(save_folder, "mean_corr_calcium_ipsi", ".pdf", sep = ""))

  plot_line_corrs_across_age(data = dev,y = "mean_corr_spikes_contra", ylab = "mean corr Coeff. (spikes)")
  ggsave(paste(save_folder, "mean_corr_spikes_contra", ".pdf", sep = ""))

  plot_line_corrs_across_age(data = dev,y = "mean_corr_calcium_contra", ylab = "mean corr Coeff. (cal)")
  ggsave(paste(save_folder, "mean_corr_calcium_contra", ".pdf", sep = ""))
}

pipline <- function(deconvolution_method = "BCL", remove_intro = "False", results = "../BCL_results/Correlation_results/BCL_intro_not_removed/"){
  o_grav_7 <- corrs_all_fish_in_group("WT_GR_7_dpf/",genotype = "WT", age = "7_dpf", rearing_condition = "GR", deconvolution_method = deconvolution_method , remove_intro = remove_intro)
  o_norm_7 <- corrs_all_fish_in_group("WT_NR_7_dpf/", genotype = "WT", age = "7_dpf", rearing_condition = "NR", deconvolution_method = deconvolution_method, remove_intro = remove_intro)

  o_grav_5 <- corrs_all_fish_in_group("WT_GR_5_dpf/", genotype = "WT", age = "5_dpf", rearing_condition = "GR", deconvolution_method = deconvolution_method, remove_intro = remove_intro)
  o_norm_5 <- corrs_all_fish_in_group("WT_NR_5_dpf/", genotype = "WT", age = "5_dpf", rearing_condition = "NR", deconvolution_method = deconvolution_method, remove_intro = remove_intro)

  o_grav_3 <- corrs_all_fish_in_group("WT_GR_3_dpf/", genotype = "WT", age = "3_dpf", rearing_condition = "GR", deconvolution_method = deconvolution_method, remove_intro = remove_intro)
  o_norm_3 <- corrs_all_fish_in_group("WT_NR_3_dpf/", genotype = "WT", age = "3_dpf", rearing_condition = "NR", deconvolution_method = deconvolution_method, remove_intro = remove_intro)

plot_correltaions_by_age(save_folder = results)
plot_corrs_across_age(save_folder =  results)
}

pipline(results = "../BCL_results/Correlation_results/BCL_intro_not_removed/")
pipline(deconvolution_method = "BCL", remove_intro = "True", results = "../BCL_results/Correlation_results/BCL_intro_removed/")