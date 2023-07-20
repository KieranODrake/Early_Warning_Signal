#' Takes data frames created on HPC using . These dataframes represent the 
#' generation of early warning signals (EWS) and assessment [of True or False Positives]
#' for different threshold levels. There are three types of data frame and a number of
#' files for each data frame relating to different non-TFPS leading indicators:
#'  - CoMix
#'  - Ct_p2_mean and median
#'  - Google mobility
#'  - PCR positivity rates
#'  - hospitalisation tests
#' 
#' This code takes all of the .rds files in a folder and combines them into three
#' dataframes depending on the names of the file

################### Create list of .rds files
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/EWS_df_non_TFPS/" )
file_list = list.files() ; file_list = subset( file_list , file_list != "desktop.ini" & file_list != "test" & file_list != "lgr_th_NA" )  
file_list_ews_results = as.data.frame( subset( file_list , grepl( "ews_results" , file_list ) ) )
file_list_wave_results = as.data.frame( subset( file_list , grepl( "wave_results_df" , file_list ) ) )
file_list_wave_results_analysis = as.data.frame( subset( file_list , grepl( "wave_results_analysis_df" , file_list ) ) )

library(data.table)

#' Merge ews_results_df_XXX.rds
ews_results_df_list = lapply( file_list_ews_results[[1]] , readRDS )
ews_results_df = data.table::rbindlist( ews_results_df_list )
#' Merge wave_results_df_XXX.rds
wave_results_df_list = lapply( file_list_wave_results[[1]] , readRDS )
wave_results_df = data.table::rbindlist( wave_results_df_list )
#' Merge wave_results_analysis_results_df_XXX.rds
wave_results_analysis_df_list = lapply( file_list_wave_results_analysis[[1]] , readRDS )
wave_results_analysis_df = data.table::rbindlist( wave_results_analysis_df_list )

#' Save files
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/EWS_df_combined_non_TFPS/" )
saveRDS( ews_results_df , "ews_results_df.rds")
rm( file_list_ews_results , ews_results_df_list )
saveRDS( wave_results_df , "wave_results_df.rds")
rm( wave_results_df_list , file_list_wave_results )
saveRDS( wave_results_analysis_df , "wave_results_analysis_df.rds")
rm( wave_results_analysis_df_list , file_list_wave_results_analysis )
gc()