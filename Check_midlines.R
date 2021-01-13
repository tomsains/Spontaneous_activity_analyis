# Title     : TODO
# Objective : TODO
# Created by: thomas_sainsbury
# Created on: 27/12/2020


library(manipulate)
library(stringr)

list_all_folders <- function(folder) {
  return(list.files(pattern = "WT*"))
}

list_all_fish_in_dir <- function(folder ){
    list_fish = list.files(path = folder, pattern = "*_all_cells_centers_cut_ordered.dat")
    #print(list_fish)
    return(list_fish)
}

find_midline <- function(centers){
  plot_points <- function(data, a, b){
    smoothScatter(data[,1],data[,2], transformation = function(x) x^.25)
    points(data[,1],data[,2], cex= 0.2)
    abline(a = a, b = b, col ="red", lwd =5)
  }
  manipulate(plot_points(centers [,1:2], a, b), a =slider(-100,500), b = slider(-2,-0.2))
}

find_all_midlines_in_folder <- function(folder){
  list_fish <- list_all_fish_in_dir(folder = folder)
  for (i in 1:length(list_fish)) {
    print(list_fish [[i]])
    centers <- as.matrix(fread(paste(folder, list_fish [[i]], sep ="")))
    find_midline(centers = centers)
    line <- readline()
  }
}



find_all_midlines_in_folder("WT_GR_3_dpf/")
write.table(x = data.frame(data_set_name = c("180810_grav_3dpf_h2b_GCAMP6_sa_f2_00003_scaled_aligned", "190304_grav_3dpf_h2b_GCAMP6_sa_F1_00001_scaled_aligned","190304_grav_3dpf_h2b_GCAMP6_sa_F3_00002_scaled_aligned", "190513_grav_3dpf_h2b_GCAMP6_sa_F2_00001_scaled_aligned"),
           intercept = c(300,280,186, 310),  slope = c(-1.26, -1.08, -0.43, -1.28)),file = "WT_GR_3_dpf/midlines.dat")






find_all_midlines_in_folder("WT_GR_5_dpf/")
write.table(x = data.frame(data_set_name = c("190306_grav_5dpf_h2b_GCAMP6_sa_F1_00002_scaled_aligned", "190306_grav_5dpf_h2b_GCAMP6_sa_F2_00002_scaled_aligned",
                                             "190312_grav_5dpf_h2b_GCAMP6_sa_F1_00002_scaled_aligned", "190312_grav_5dpf_h2b_GCAMP6_sa_F2_00001_scaled_aligned",
                                             "190515_grav_5dpf_h2b_GCAMP6_sa_F1_00001_scaled_aligned"),
                           intercept = c(220,197,399,255, 205),  slope = c(-0.91, -0.69,-1.96, -0.77, -0.73)),file = "WT_GR_5_dpf/midlines.dat")






find_all_midlines_in_folder("WT_GR_7_dpf/")
write.table(x = data.frame(data_set_name = c("180530_WT_grav_h2b_gc6s_7dpf_f1_sa__00001_scaled_aligned", "180530_WT_grav_h2b_gc6s_7dpf_f2_sa__00001_scaled_aligned",
                                             "180607_WT_grav_h2b_gc6s_7dpf_f1_sa__00001_scaled_aligned", "180607_WT_grav_h2b_gc6s_7dpf_f2_sa__00002_scaled_aligned",
                                             "180622_WT_grav_h2b_gc6s_7dpf_f2_sa__00006_scaled_aligned", "180724_WT_grav_h2b_gc6s_7dpf_f1_sa__00001_scaled_aligned",
                                             "180724_WT_grav_h2b_gc6s_7dpf_f3_sa__00001_scaled_aligned"),
                           intercept = c(318,313,263,230,399,369,407),  slope = c(-1.41, -1.29,-1.05, -0.73, -2,-1.73,-1.97)),file = "WT_GR_7_dpf/midlines.dat")

find_all_midlines_in_folder("WT_NR_7_dpf/")
write.table(x = data.frame(data_set_name = c("180220_WT_h2b_gc6s_7dpf_f1_sa_aligned", "180220_WT_h2b_gc6s_7dpf_f2_sa_aligned",
                                             "180220_WT_h2b_gc6s_7dpf_f3_sa_aligned", "180323_WT_h2b_gc6s_7dpf_f1_sa_aligned",
                                             "180503_WT_h2b_gc6s_7dpf_f1_sa_aligned", "180510_WT_h2b_gc6s_7dpf_f1_sa_aligned",
                                             "180510_WT_h2b_gc6s_7dpf_f2_sa_aligned", "180510_WT_h2b_gc6s_7dpf_f3_sa_aligned"),
                           intercept = c(253,313,238,334,387,344,240,356),  slope = c(-1.0,-1.2,-0.85,-1.51,-1.85,-1.48,-0.99,-1.6)),file = "WT_NR_7_dpf/midlines.dat")



find_all_midlines_in_folder("WT_NR_5_dpf/")
write.table(x = data.frame(data_set_name = c("180129_WT_h2b_gc6s_5dpf_f2_sa_aligned", "180129_WT_h2b_gc6s_5dpf_f3_sa_aligned",
                                             "180321_WT_h2b_gc6s_5dpf_f1_sa_aligned", "180321_WT_h2b_gc6s_5dpf_f2_sa_aligned",
                                             "180328_WT_h2b_gc6s_5dpf_f2_sa_aligned", "180717_WT_h2b_gc6s_5dpf_f1_sa_aligned",
                                             "180717_WT_h2b_gc6s_5dpf_f3_sa_aligned","180718_WT_h2b_gc6s_5dpf_f1_sa_aligned"),
                           intercept = c(374,392,424,288,326,218,315,417),  slope = c(-1.8,-2,-2,-1.15,-1.48,-0.74,-1.51,-2)),file = "WT_NR_5_dpf/midlines.dat")



NR_3 <- find_all_midlines_in_folder("WT_NR_3_dpf/")
write.table(x = data.frame(data_set_name = c("180209_WT_h2b_gc6s_3dpf_f1_sa_aligned","180219_WT_h2b_gc6s_3dpf_f2_sa_aligned",
                                             "180219_WT_h2b_gc6s_3dpf_f4_sa_aligned", "180319_WT_h2b_gc6s_3dpf_f2_sa_aligned",
                                             "180625_WT_h2b_gc6s_3dpf_f1_sa_aligned", "180625_WT_h2b_gc6s_3dpf_f3_sa_aligned"),
                           intercept = c(402,293,324,281,271,250),  slope = c(-1.9,-1.23,-1.38,-1.136,-1.2,-1.03)),file = "WT_NR_3_dpf/midlines.dat")





find_all_midlines_in_folder("WT_GR_5_dpf/")
write.table(x = data.frame(data_set_name = c("190306_grav_5dpf_h2b_GCAMP6_sa_F1_00002_scaled_aligned", "190306_grav_5dpf_h2b_GCAMP6_sa_F2_00002_scaled_aligned",
                                             "190312_grav_5dpf_h2b_GCAMP6_sa_F1_00002_scaled_aligned", "190312_grav_5dpf_h2b_GCAMP6_sa_F2_00001_scaled_aligned",
                                             "190515_grav_5dpf_h2b_GCAMP6_sa_F1_00001_scaled_aligned"),
                           intercept = c(220,197,399,255, 205),  slope = c(-0.91, -0.69,-1.96, -0.77, -0.73)),file = "WT_GR_5_dpf/midlines.dat")


