#' Splitting the counting of false positive EWS between those generated before 
#' the earliest true psoitive and those generated after the earliest true 
#' positive for each parameter set
#' 

#' array number from High Performance Computing (HPC) array job number
params_run <- commandArgs( trailingOnly = TRUE )
#params_run="2 1 100000"
params <- strsplit( params_run , "," )
wave_n <- as.numeric( params[[1]][1] )
row_start <- as.numeric( params[[1]][2] )
row_end <- as.numeric( params[[1]][3] )

message( "Wave number: ", wave_n , ". " , "Row start: " , row_start , ". " , "Row end: " , row_end , "." )

library(data.table)
library(dplyr)

#' Load environment pre-prepared using '1_false_positive_split_prep_for_HPC.R'
#setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/false_positive_split")
load("false_positive_split_env_HPC_start_v2.RData" )

#' Run function for each wave and save output
ptm <- proc.time()
output_list = fp_split( wave_results = wave_results_dt
                      , earliest_tp = earliest_tp_dt[ row_start : row_end , ]
                      , wave_number = wave_n
                      )
proc.time() - ptm

filename_list <- paste0("output_list_w",wave_n,"_",row_start,"_",row_end,".rds")
saveRDS( output_list , filename_list )
output_df = dplyr::bind_rows(output_list)
filename_df <- paste0("output_df_w",wave_n,"_",row_start,"_",row_end,".rds")
saveRDS( output_df , filename_df )
