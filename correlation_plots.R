
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

load_corr <- function(folder = "WT_GR_7_dpf/correlations/", data_set = "180530_WT_grav_h2b_gc6s_7dpf_f1_sa__00001_scaled_aligned", deconvolution_method = "BCL", remove_intro = "True", iter = 100) {
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
  
  same_side_cor_spikes <- by(corrs$spikes_corr [(corrs$spikes_corr > corrs$null_corrs) & (corrs$same_side == 1)], cut(corrs$distance_between_cells[(corrs$spikes_corr  > corrs$null_corrs) & (corrs$same_side == 1)]*1.3, breaks =seq(0, 300, 50)), mean, na.rm = TRUE)
  same_side_cor_cal <- by(corrs$cal_corr [(corrs$spikes_corr  > corrs$null_corrs) & (corrs$same_side == 1)], cut(corrs$distance_between_cells[(corrs$spikes_corr  > corrs$null_corrs) & (corrs$same_side == 1)]*1.3, breaks =seq(0, 300, 50)), mean, na.rm = TRUE)
  
  same_side_cor_spikes_no_thresh <- by(corrs$spikes_corr [(corrs$spikes_corr  > 0) &(corrs$same_side == 1)], cut(corrs$distance_between_cells[(corrs$spikes_corr  > 0) & (corrs$same_side == 1)]*1.3, breaks =seq(0, 300, 50)), mean, na.rm = TRUE)
  same_side_cor_cal_no_thresh <- by(corrs$cal_corr [(corrs$spikes_corr  > 0) & (corrs$same_side == 1)], cut(corrs$distance_between_cells[(corrs$spikes_corr > 0) & (corrs$same_side == 1)]*1.3, breaks =seq(0, 300, 50)), mean, na.rm = TRUE)
  return(c(as.vector(same_side_cor_spikes), as.vector(same_side_cor_cal), as.vector(same_side_cor_spikes_no_thresh), as.vector(same_side_cor_cal_no_thresh)))
}

cut_corrs_by_dist_difference <- function(corrs = corrs){
  same_side_cor_spikes <- by((corrs$spikes_corr  - corrs$null_corrs) [(corrs$spikes_corr > corrs$null_corrs) & (corrs$same_side == 1)], cut(corrs$distance_between_cells[(corrs$spikes_corr  > corrs$null_corrs) & (corrs$same_side == 1)]*1.3, breaks =seq(0, 300, 50)), mean, na.rm = TRUE)
  same_side_cor_cal <- by((corrs$spikes_corr  - corrs$null_corrs) [(corrs$spikes_corr  > corrs$null_corrs) & (corrs$same_side == 1)], cut(corrs$distance_between_cells[(corrs$spikes_corr  > corrs$null_corrs) & (corrs$same_side == 1)]*1.3, breaks =seq(0, 300, 50)), mean, na.rm = TRUE)
  
  return(c(as.vector(same_side_cor_spikes), as.vector(same_side_cor_cal)))
}


get_lower_tri_mat <- function(mat){
  return(mat[lower.tri(mat)])
}


pairwise_correlation_summaries <- function(folder = "WT_GR_7_dpf/correlations/", data_set = "180530_WT_grav_h2b_gc6s_7dpf_f1_sa__00001_scaled_aligned", deconvolution_method, remove_intro){
  corrs <- load_corr(folder  = folder, data_set = data_set, remove_intro = remove_intro, deconvolution_method = deconvolution_method)
  data <- get_lowers_tri_for_corrs(corrs)
  dat <- as.data.frame(cbind(mean(data$spikes_corr [(data$spikes_corr > data$null_corrs) & (data$same_side == 1)], na.rm = TRUE),
               mean(data$cal_corr [(data$spikes_corr > data$null_corrs) & (data$same_side == 1)], na.rm = TRUE),
               mean(data$spikes_corr [(data$spikes_corr > data$null_corrs) & (data$same_side == -1)], na.rm = TRUE),
               mean(data$cal_corr [(data$spikes_corr > data$null_corrs) & (data$same_side == -1)], na.rm = TRUE),
               median(data$spikes_corr [(data$spikes_corr > data$null_corrs) & (data$same_side == 1)], na.rm = TRUE),
               median(data$cal_corr [(data$spikes_corr > data$null_corrs) & (data$same_side == 1)], na.rm = TRUE),
               median(data$spikes_corr [(data$spikes_corr > data$null_corrs) & (data$same_side == -1)], na.rm = TRUE),
               median(data$cal_corr [(data$spikes_corr > data$null_corrs) & (data$same_side == -1)], na.rm = TRUE),
               mean(data$spikes_corr [(data$spikes_corr > 0) & (data$same_side == 1)], na.rm = TRUE),
               mean(data$cal_corr [(data$spikes_corr > 0) & (data$same_side == 1)], na.rm = TRUE),
               mean(data$spikes_corr [(data$spikes_corr > 0) & (data$same_side == -1)], na.rm = TRUE),
               mean(data$cal_corr [(data$spikes_corr > 0) & (data$same_side == -1)], na.rm = TRUE),
               median(data$spikes_corr [(data$spikes_corr > 0) & (data$same_side == 1)], na.rm = TRUE),
               median(data$cal_corr [(data$spikes_corr > 0) & (data$same_side == 1)], na.rm = TRUE),
               median(data$spikes_corr [(data$spikes_corr > 0) & (data$same_side == -1)], na.rm = TRUE),
               median(data$cal_corr [(data$spikes_corr > 0) & (data$same_side == -1)], na.rm = TRUE),
               t(cut_corrs_by_dist(data)), t(cut_corrs_by_dist_difference(data))))
  colnames(dat) <- c("mean_corr_spikes_ipsi", "mean_corr_calcium_ipsi", "mean_corr_spikes_contra", "mean_corr_calcium_contra",
                     "median_corr_spikes_ipsi", "median_corr_calcium_ipsi", "median_corr_spikes_conta", "median_corr_spikes_contra",
                     "NT_mean_corr_spikes_ipsi", "NT_mean_corr_calcium_ipsi", "NT_mean_corr_spikes_contra", "NT_mean_corr_calcium_contra",
                     "NT_median_corr_spikes_ipsi", "NT_median_corr_calcium_ipsi", "NT_median_corr_spikes_conta", "NT_median_corr_spikes_contra",
                     sprintf("Distance_spike_corr_%03d", seq(0,250,50)),sprintf("Distance_cal_corr_%03d", seq(0,250,50)),sprintf("NT_Distance_spike_corr_%03d", seq(0,250,50)),sprintf("NT_Distance_cal_corr_%03d", seq(0,250,50)),
                     sprintf("diff_Distance_spike_corr_%03d", seq(0,250,50)),sprintf("diff_Distance_cal_corr_%03d", seq(0,250,50)))
 # print(dat)
  
  
  dat$ipsi_contra_ratio <- dat$mean_corr_spikes_ipsi/dat$mean_corr_spikes_contra
  #print("worked")
  dat$mean_corr_spikes_distance_difference <- mean(c(dat$Distance_spike_corr_000, dat$Distance_spike_corr_050))/mean(c(dat$Distance_spike_corr_250 + dat$Distance_spike_corr_200))
  
  
  dat$NT_ipsi_contra_ratio <- dat$NT_mean_corr_spikes_ipsi/dat$NT_mean_corr_spikes_contra
  #print("worked")
  dat$NT_mean_corr_spikes_distance_difference <- mean(c(dat$NT_Distance_spike_corr_000, dat$NT_Distance_spike_corr_050))/mean(c(dat$NT_Distance_spike_corr_250 + dat$NT_Distance_spike_corr_200))
  #dat$NT_mean_spike_difference <- c(((dat$NT_Distance_spike_corr_000 + dat$NT_Distance_spike_corr_100 + dat$NTDistance_spike_corr_050)/3)/((dat$NT_Distance_spike_corr_250 + dat$NT_Distance_spike_corr_200)/2))
  # 
  dat$short_range_correlations <- mean(c(dat$Distance_spike_corr_000, dat$Distance_spike_corr_050))
  dat$long_range_correlations <- mean(c(dat$Distance_spike_corr_250, dat$Distance_spike_corr_200))
  dat$long_short_ratio_spikes <- dat$short_range_correlations/dat$long_range_correlations
  # 
  dat$NT_short_range_correlations <- mean(c(dat$NT_Distance_spike_corr_000, dat$NTDistance_spike_corr_050))
  dat$NT_long_range_correlations <- mean(c(dat$NT_Distance_spike_corr_250, dat$NT_Distance_spike_corr_200))
  dat$NT_long_short_ratio_spikes <- dat$NT_short_range_correlations/dat$NT_long_range_correlations
  # # 
  dat$NT_short_range_correlations_cal <- mean(c(dat$NT_Distance_cal_corr_000, dat$NT_Distance_cal_corr_050))
  dat$NT_long_range_correlations_cal <- mean(c(dat$NT_Distance_cal_corr_200, dat$NT_Distance_cal_corr_250))
  dat$NT_long_short_ratio_cal <- dat$NT_short_range_correlations_cal/dat$NT_long_range_correlations_cal
  # # 
  dat$short_range_correlations_cal <- mean(c(dat$Distance_cal_corr_000, dat$Distance_cal_corr_050))
  dat$long_range_correlations_cal <- mean(c(dat$Distance_cal_corr_200, dat$Distance_cal_corr_250))
  dat$long_short_ratio_cal <- dat$short_range_correlations_cal/dat$long_range_correlations_cal
  
  
  
  
  
  
 # dat$diff_ipsi_contra_ratio <- dat$diff_mean_corr_spikes_ipsi/dat$diff_mean_corr_spikes_contra
  #print("worked")
  dat$diff_mean_corr_spikes_distance_difference <- mean(c(dat$diff_Distance_spike_corr_000, dat$diff_Distance_spike_corr_050))/mean(c(dat$diff_Distance_spike_corr_250 + dat$diff_Distance_spike_corr_200))
  #dat$diff_mean_spike_difference <- c(((dat$diff_Distance_spike_corr_000 + dat$diff_Distance_spike_corr_100 + dat$NTDistance_spike_corr_050)/3)/((dat$diff_Distance_spike_corr_250 + dat$diff_Distance_spike_corr_200)/2))
  # 

  # 
  dat$diff_short_range_correlations <- mean(c(dat$diff_Distance_spike_corr_000, dat$NTDistance_spike_corr_050))
  dat$diff_long_range_correlations <- mean(c(dat$diff_Distance_spike_corr_250, dat$diff_Distance_spike_corr_200))
  dat$diff_long_short_ratio_spikes <- dat$diff_short_range_correlations/dat$diff_long_range_correlations
  # # 

  # # 
  dat$short_range_correlations_cal <- mean(c(dat$Distance_cal_corr_000, dat$Distance_cal_corr_050))
  dat$long_range_correlations_cal <- mean(c(dat$Distance_cal_corr_200, dat$Distance_cal_corr_250))
  dat$long_short_ratio_cal <- dat$short_range_correlations_cal/dat$long_range_correlations_cal
  # 
  return(dat)
  #print("returned")
}

list_all_fish_in_dir <- function(folder, deconvolution_method = "BCL", remove_intro ){
  list_fish = list.files(path = paste(folder, sep =""), pattern = paste(deconvolution_method, "_cell_side_removed_intro-", remove_intro, ".dat", sep =""))
  print(list_fish)
  return(lapply(list_fish, str_replace, pattern = paste(deconvolution_method, "_cell_side_removed_intro-", remove_intro, ".dat", sep =""), replacement = ""))
}

corrs_all_fish_in_group <- function(folder, deconvolution_method = "BCL", remove_intro = "True"){
  folder <- paste(folder, "correlations/", sep ="")
  data_sets <- list_all_fish_in_dir(folder=folder, deconvolution_method = deconvolution_method, remove_intro = remove_intro)
  
  cor_for_group <- vector("list", length(data_sets))
  for (i in 1:length(data_sets)){
    print(data_sets [[i]])
    cor_for_group [[i]] <- pairwise_correlation_summaries(folder = folder, data_set = data_sets[i],deconvolution_method=deconvolution_method, remove_intro = remove_intro)
  }
  cor_for_group <- as.data.frame(do.call(rbind, cor_for_group))
  group_stats <- data.frame("data_set_names" = unlist(data_sets), "Genotype" = rep(substring(folder, 1, 2), length(data_sets)), "Age" = rep(substring(folder, 7, 7), length(data_sets)), "Rearing_condition" = rep(substring(folder, 4, 5), length(data_sets)))
  dat <-cbind(group_stats,cor_for_group)
  return(dat)
}




plot_correlations_for_each_age <- function(dev, age = "5_dpf", save_folder = "../BCL_results/Correlation_results/BCL_intro_removed/") {
  
  if (age == "3_dpf") {
    combine <- dev [dev$Age  == 3,]
  }
  if (age == "5_dpf") {
    combine <- dev [dev$Age  == 5,]
  }
  if (age == "7_dpf") {
    combine <- dev [dev$Age  == 7,]
  }
  
  combine %>% select(Rearing_condition, sprintf("Distance_spike_corr_%03d", seq(0,250,50))) %>%
    pivot_longer(., cols = c(sprintf("Distance_spike_corr_%03d", seq(0,250,50))), names_to = "Var", values_to = "Val") %>%
    ggplot(aes(x = Var, y = Val, fill = Rearing_condition)) +
    geom_boxplot() + xlab("Distance (um)") + ylab("Correlation Coef.") + scale_x_discrete(labels=c("0-50", "50-100", "100-150", "150-200", "200-250", "250-300"))
  
  ggsave(paste(save_folder, "Distance_spikes_corrs_", age, ".pdf", sep = ""))
  
  combine %>% select(Rearing_condition, sprintf("diff_Distance_spike_corr_%03d", seq(0,250,50))) %>%
    pivot_longer(., cols = c(sprintf("diff_Distance_spike_corr_%03d", seq(0,250,50))), names_to = "Var", values_to = "Val") %>%
    ggplot(aes(x = Var, y = Val, fill = Rearing_condition)) +
    geom_boxplot() + xlab("Distance (um)") + ylab("Correlation Coef.") + scale_x_discrete(labels=c("0-50", "50-100", "100-150", "150-200", "200-250", "250-300"))
  
  ggsave(paste(save_folder, "diff_Distance_spikes_corrs_", age, ".pdf", sep = ""))
  
  
  combine %>% select(Rearing_condition, sprintf("NT_Distance_spike_corr_%03d", seq(0,250,50))) %>%
    pivot_longer(., cols = c(sprintf("NT_Distance_spike_corr_%03d", seq(0,250,50))), names_to = "Var", values_to = "Val") %>%
    ggplot(aes(x = Var, y = Val, fill = Rearing_condition)) +
    geom_boxplot() + xlab("Distance (um)") + ylab("Correlation Coef.") + scale_x_discrete(labels=c("0-50", "50-100", "100-150", "150-200", "200-250", "250-300"))
  
  ggsave(paste(save_folder, "NT_Distance_spikes_corrs_", age, ".pdf", sep = ""))
  
  combine %>% select(Rearing_condition, sprintf("Distance_cal_corr_%03d", seq(0,250,50))) %>%
    pivot_longer(., cols = c(sprintf("Distance_cal_corr_%03d", seq(0,250,50))), names_to = "Var", values_to = "Val") %>%
    ggplot(aes(x = Var, y = Val, fill = Rearing_condition)) +
    geom_boxplot() + xlab("Distance (um)") + ylab("Correlation Coef.") + scale_x_discrete(labels=c("0-50", "50-100", "100-150", "150-200", "200-250", "250-300"))
  ggsave(paste(save_folder, "Distance_cal_corrs_", age, ".pdf", sep = ""))
  
  combine %>% select(Rearing_condition, sprintf("NT_Distance_cal_corr_%03d", seq(0,250,50))) %>%
    pivot_longer(., cols = c(sprintf("NT_Distance_cal_corr_%03d", seq(0,250,50))), names_to = "Var", values_to = "Val") %>%
    ggplot(aes(x = Var, y = Val, fill = Rearing_condition)) +
    geom_boxplot() + xlab("Distance (um)") + ylab("Correlation Coef.") + scale_x_discrete(labels=c("0-50", "50-100", "100-150", "150-200", "200-250", "250-300"))
  ggsave(paste(save_folder, "NT_Distance_cal_corrs_", age, ".pdf", sep = ""))
  
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

plot_correltaions_by_age <- function(dev, Save_folder){
  plot_correlations_for_each_age(dev = dev, age = "3_dpf", save_folder = Save_folder)
  plot_correlations_for_each_age(dev = dev, age = "5_dpf", save_folder  = Save_folder)
  plot_correlations_for_each_age(dev = dev, age = "7_dpf", save_folder = Save_folder)
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
  data %>% select("Age", "Rearing_condition", all_of(y)) %>% ggplot(aes_string(x = "Age", y= y, group = "Rearing_condition", col = "Rearing_condition"))  + ylab(ylab) +stat_summary(fun=mean, geom="line", lwd = 1)  +
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.2,lwd = 1) +
    stat_summary(fun=mean, geom="point", size = 5) + xlab("Age (dpf)") + theme_classic(base_size = 40) +
    scale_color_manual(values = col) + theme(legend.position="none") + geom_jitter(width = 0.05, alpha =0.5, size = 3)
}

plot_corrs_across_age <- function(dev, save_folder = "../BCL_results/Correlation_results/BCL_intro_removed/"){
  #dev <- as.data.frame(rbind(o_norm_3, o_grav_3, o_norm_5, o_norm_5, o_grav_5, o_norm_7, o_grav_7))
  
  plot_line_corrs_across_age(data = dev, y = "short_range_correlations")
  ggsave(paste(save_folder, "Short_range_corrs_across_development", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "long_range_correlations", ylab = "Long range corr coeff")
  ggsave(paste(save_folder, "Long_range_corrs_across_development", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "long_short_ratio_spikes", ylab = "short/long range corrs")
  ggsave(paste(save_folder, "Long_short_corrs_ratio_across_development", ".pdf", sep = ""))
  
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
  
  
  # No threshold on circular permuted corrs applied
  
  plot_line_corrs_across_age(data = dev, y = "NT_short_range_correlations")
  ggsave(paste(save_folder, "NT_Short_range_corrs_across_development", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "NT_long_range_correlations", ylab = "Long range corr coeff")
  ggsave(paste(save_folder, "NT_Long_range_corrs_across_development", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "NT_long_short_ratio_spikes", ylab = "short/long range corrs")
  ggsave(paste(save_folder, "NT_Long_short_corrs_ratio_across_development", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev, y = "NT_short_range_correlations_cal")
  ggsave(paste(save_folder, "NT_Short_range_corrs_across_development_cal", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "NT_long_range_correlations_cal", ylab = "Long range corr coeff")
  ggsave(paste(save_folder, "NT_Long_range_corrs_across_development_cal", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "NT_long_short_ratio_cal", ylab = "short/long range corrs")
  ggsave(paste(save_folder, "Long_short_corrs_ratio_across_development_cal", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "NT_mean_corr_spikes_ipsi", ylab = "mean corr Coeff. (spikes)")
  ggsave(paste(save_folder, "mean_corr_spikes_ipsi", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "NT_mean_corr_calcium_ipsi", ylab = "mean corr Coeff. (cal)")
  ggsave(paste(save_folder, "NT_mean_corr_calcium_ipsi", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "NT_mean_corr_spikes_contra", ylab = "mean corr Coeff. (spikes)")
  ggsave(paste(save_folder, "NT_mean_corr_spikes_contra", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "NT_mean_corr_calcium_contra", ylab = "mean corr Coeff. (cal)")
  ggsave(paste(save_folder, "NT_mean_corr_calcium_contra", ".pdf", sep = ""))
  
  
  # diffs_corrs
  plot_line_corrs_across_age(data = dev, y = "diff_short_range_correlations")
  ggsave(paste(save_folder, "diff_Short_range_corrs_across_development", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "diff_long_range_correlations", ylab = "Long range corr coeff")
  ggsave(paste(save_folder, "diff_Long_range_corrs_across_development", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "diff_long_short_ratio_spikes", ylab = "short/long range corrs")
  ggsave(paste(save_folder, "diff_Long_short_corrs_ratio_across_development", ".pdf", sep = ""))
  
 
  plot_line_corrs_across_age(data = dev,y = "diff_long_range_correlations_cal", ylab = "Long range corr coeff")
  ggsave(paste(save_folder, "diff_Long_range_corrs_across_development_cal", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "diff_long_short_ratio_cal", ylab = "short/long range corrs")
  ggsave(paste(save_folder, "Long_short_corrs_ratio_across_development_cal", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "diff_mean_corr_spikes_ipsi", ylab = "mean corr Coeff. (spikes)")
  ggsave(paste(save_folder, "mean_corr_spikes_ipsi", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "diff_mean_corr_calcium_ipsi", ylab = "mean corr Coeff. (cal)")
  ggsave(paste(save_folder, "diff_mean_corr_calcium_ipsi", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "diff_mean_corr_spikes_contra", ylab = "mean corr Coeff. (spikes)")
  ggsave(paste(save_folder, "diff_mean_corr_spikes_contra", ".pdf", sep = ""))
  
  plot_line_corrs_across_age(data = dev,y = "diff_mean_corr_calcium_contra", ylab = "mean corr Coeff. (cal)")
  ggsave(paste(save_folder, "diff_mean_corr_calcium_contra", ".pdf", sep = ""))
}

pipline <- function(deconvolution_method = "BCL", remove_intro = "True", results = "../BCL_results/Correlation_results/BCL_intro_not_removed/"){
  o_grav_7 <- corrs_all_fish_in_group("WT_GR_7_dpf/", deconvolution_method = deconvolution_method , remove_intro = remove_intro)
  o_norm_7 <- corrs_all_fish_in_group("WT_NR_7_dpf/", deconvolution_method = deconvolution_method, remove_intro = remove_intro)
  
  o_grav_5 <- corrs_all_fish_in_group("WT_GR_5_dpf/", deconvolution_method = deconvolution_method, remove_intro = remove_intro)
  o_norm_5 <- corrs_all_fish_in_group("WT_NR_5_dpf/", deconvolution_method = deconvolution_method, remove_intro = remove_intro)
  
  o_grav_3 <- corrs_all_fish_in_group("WT_GR_3_dpf/",deconvolution_method = deconvolution_method, remove_intro = remove_intro)
  o_norm_3 <- corrs_all_fish_in_group("WT_NR_3_dpf/", deconvolution_method = deconvolution_method, remove_intro = remove_intro)
  dev <- as.data.frame(rbind(o_norm_3, o_grav_3, o_norm_5, o_grav_5, o_norm_7, o_grav_7))
  plot_correltaions_by_age(dev = dev,Save_folder = results)
  plot_corrs_across_age(dev= dev, save_folder = results)
  return(dev)
}



#pipline(deconvolution_method = "BCL", remove_intro = "True", results = "../BCL_results/Correlation_results/BCL_intro_removed/")

#pipline(deconvolution_method = "BCL", remove_intro = "False", results = "../BCL_results/Correlation_results/BCL_intro_not_removed/")


setwd("/media/thomas_sainsbury/spont//Spontaneous_activity_experiments/OASIS_output/")
pipline(deconvolution_method = "AR1", remove_intro = "True", results = "../OASIS_results/Correlation_results/AR1_intro_removed/")
pipline(deconvolution_method = "AR1", remove_intro = "False", results = "../OASIS_results/Correlation_results/AR1_intro_not_removed/")


pipline(deconvolution_method = "estimated", remove_intro = "True", results = "../OASIS_results/Correlation_results/estimated_intro_removed/")
pipline(deconvolution_method = "estimated", remove_intro = "False", results = "../OASIS_results/Correlation_results/estimated_intro_not_removed/")
