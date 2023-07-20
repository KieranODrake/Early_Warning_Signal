#' Merge false positive split dataframes produced using HPC
#library(dplyr)
library(data.table)
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/false_positive_split_2023_03/output" )
file_list = list.files() ; file_list = subset( file_list , file_list != "desktop.ini" ) 

df_list = list()
for (i in 2:8){
  #' Create list of dataframes to merge for each wave
  temp = file_list[ grep( paste0("output_df_w",i) , file_list ) ]
  #print(temp)
  new_name = paste0( "file_list_w",i )
  assign( new_name , temp )

  #' Merge dataframes for each wave
  fp_split_df_temp = data.frame()
  for (j in 1 : length( temp ) ){
    temp_df <- readRDS( temp[ j ] )  #temp_df <- readRDS( paste0(folder,file_list[ i ] ) )
    fp_split_df_temp <- rbind( fp_split_df_temp , temp_df )
  }
  new_name_df = paste0( "fp_split_df_w" , i )
  assign( new_name_df , fp_split_df_temp )
  df_list = append( df_list , new_name_df) #eval( as.name( new_name_df ) ) )
}

#' Then merge combined dataframes for each wave
library(tidyverse)
#merge all data frames together
#fp_split_df = eval( as.name( df_list ) ) %>% reduce( full_join , by='variables_code')
df_list <- list( fp_split_df_w2 , fp_split_df_w3 , fp_split_df_w4 , fp_split_df_w5 , fp_split_df_w6 , fp_split_df_w7 , fp_split_df_w8 )
fp_split_df =  df_list %>% reduce( full_join , by='variables_code')

rm(df_list, fp_split_df_temp, temp_df)
rm(fp_split_df_w2,fp_split_df_w3,fp_split_df_w4,fp_split_df_w5,fp_split_df_w6,fp_split_df_w7,fp_split_df_w8)
gc()
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/false_positive_split_2023_03/" )
saveRDS( fp_split_df , "fp_split_df.rds" )

#' Then merge into analysis for paper
#' in wave_results_analysis_Feb_2023