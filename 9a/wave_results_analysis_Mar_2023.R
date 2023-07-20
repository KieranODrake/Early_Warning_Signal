#' Extract information from wave_results_analysis_df
#' Author: Kieran Drake
#' Date: November 2022
#'
#' Read in dataframe containing information on early warning signals (as calculated in EWS_calc_threshold.R) 
#' for each set of scan variables (min age, max age, min descendants) and filter variable (LGR p-value, 
#' parent/sub-cluster replacement LGR threshold) and EWS threshold
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined_2023_03_20/")
wave_results_analysis_df = readRDS( "wave_results_analysis_df.rds" )

#' wave_results_analysis_df has wave information in columns with single columns for EWS information,
#' whereas it can be searched more easily if the EWS information is combined with the wave number.
#' However, this dataframe has many duplicate sets of variables/filters (one set for each wave) and 
#' will have NAs in the new columns where the wave number doesn't match. This is handled by splitting 
#' the dataframe later and recombining to remove the NAs.

#' Initialise reshaped dataframe
var_df = data.frame( "leading_indicator_type" = rep(NA,1199651) 
                     , "TFPS_cluster_min_age" = rep(NA,1199651)
                     , "TFPS_cluster_max_age" = rep(NA,1199651)
                     , "TFPS_cluster_min_descendants" = rep(NA,1199651)
                     , "LGR_p_value_threshold" = rep(NA,1199651)
                     , "parent_sub_lgr_threshold" = rep(NA,1199651)
                     , "EWS_threshold" = rep(NA,1199651)
                     , "w2_earliest_tp_EWS_date"= rep(NA,1199651) , "w2_lead_lag_days"= rep(NA,1199651) , "w2_EWS_false_positive_n"= rep(NA,1199651)
                     , "w3_earliest_tp_EWS_date"= rep(NA,1199651) , "w3_lead_lag_days"= rep(NA,1199651) , "w3_EWS_false_positive_n"= rep(NA,1199651)
                     , "w4_earliest_tp_EWS_date"= rep(NA,1199651) , "w4_lead_lag_days"= rep(NA,1199651) , "w4_EWS_false_positive_n"= rep(NA,1199651)
                     , "w5_earliest_tp_EWS_date"= rep(NA,1199651) , "w5_lead_lag_days"= rep(NA,1199651) , "w5_EWS_false_positive_n"= rep(NA,1199651)
                     , "w6_earliest_tp_EWS_date"= rep(NA,1199651) , "w6_lead_lag_days"= rep(NA,1199651) , "w6_EWS_false_positive_n"= rep(NA,1199651)
                     , "w7_earliest_tp_EWS_date"= rep(NA,1199651) , "w7_lead_lag_days"= rep(NA,1199651) , "w7_EWS_false_positive_n"= rep(NA,1199651)
                     , "w8_earliest_tp_EWS_date"= rep(NA,1199651) , "w8_lead_lag_days"= rep(NA,1199651) , "w8_EWS_false_positive_n"= rep(NA,1199651)
)

var_df = var_df[1,]

#' Populate wave/EWS info columns
waves = na.omit( unique( wave_results_analysis_df$wave_n ) ) 
#waves = unique( wave_results_analysis_df$wave_n )
#waves = na.omit( unique( wave_results_analysis_df_inc_new_pango_stat$wave_n ) )
for ( n in waves ){
  df_filtered = subset( wave_results_analysis_df , wave_results_analysis_df$wave_n == n )
  
  message( n - 1 , " of " , length(waves), " waves complete. Rows: ", nrow(df_filtered) )
  
  #df_filtered = subset( wave_results_analysis_df_inc_new_pango_stat , wave_results_analysis_df_inc_new_pango_stat$wave_n == n )
  start_row = nrow( var_df ) + 1
  end_row = start_row + nrow( df_filtered ) - 1
  
  var_df[ start_row : end_row , 1:7 ] = df_filtered[ , c(12,7,8,9,10,11,13)]
  
  var_df[ start_row : end_row , 2+(3*n) ] = df_filtered$earliest_tp_EWS_date
  var_df[ start_row : end_row , 3+(3*n) ] = df_filtered$lead_lag_days
  var_df[ start_row : end_row , 4+(3*n) ] = df_filtered$ews_false_positive_n
}

#' Check size of var_df compared with wave_results_analysis_df
if( nrow( var_df ) != nrow( wave_results_analysis_df ) ){ message( "Dataframes have a different number of rows" ) }
'%ni%' = Negate('%in%')
View( subset( wave_results_analysis_df , wave_results_analysis_df$wave_n %ni% waves ) ) #' All of these are NA entries
unique( subset( wave_results_analysis_df , wave_results_analysis_df$wave_n %ni% waves )$wave_n )

#' Create a single value representing all of the variables/filters which can be used to filter the data
var_df$variables_code = paste0( var_df$leading_indicator_type, "_"
                                , var_df$TFPS_cluster_min_age , "_"
                                , var_df$TFPS_cluster_max_age , "_"
                                , var_df$TFPS_cluster_min_descendants , "_"
                                , var_df$LGR_p_value_threshold , "_"
                                , var_df$parent_sub_lgr_threshold , "_"
                                , var_df$EWS_threshold
)

#' Cut dataframe into one for each wave, remove NAs, then merge using the 'variables_code' and recreate variable columns
var_code_df = var_df[ , c( 29 , seq(1,7,1) )]
#' There is some repetition in the variables/filters/codes due to separate rows for each wave, so we select unique entries
var_code_df = unique( var_code_df )

#' Split main dataframe into entries by wave
w2_var_df = subset( var_df , !is.na( var_df$w2_EWS_false_positive_n ) )[ , c( 29 ,  8 ,  9 , 10 )]
w3_var_df = subset( var_df , !is.na( var_df$w3_EWS_false_positive_n ) )[ , c( 29 , 11 , 12 , 13 )]
w4_var_df = subset( var_df , !is.na( var_df$w4_EWS_false_positive_n ) )[ , c( 29 , 14 , 15 , 16 )]
w5_var_df = subset( var_df , !is.na( var_df$w5_EWS_false_positive_n ) )[ , c( 29 , 17 , 18 , 19 )]
w6_var_df = subset( var_df , !is.na( var_df$w6_EWS_false_positive_n ) )[ , c( 29 , 20 , 21 , 22 )]
w7_var_df = subset( var_df , !is.na( var_df$w7_EWS_false_positive_n ) )[ , c( 29 , 23 , 24 , 25 )]
w8_var_df = subset( var_df , !is.na( var_df$w8_EWS_false_positive_n ) )[ , c( 29 , 26 , 27 , 28 )]
#'Check if captured all (non-NA) rows
nrow( var_df )
nrow(w2_var_df) + nrow(w3_var_df) + nrow(w4_var_df) +nrow(w5_var_df) +nrow(w6_var_df) +nrow(w7_var_df) +nrow(w8_var_df) 
nrow( subset( var_df , !is.na( var_df$leading_indicator_type ) ) ) == nrow(w2_var_df) + nrow(w3_var_df) + nrow(w4_var_df) +nrow(w5_var_df) +nrow(w6_var_df) +nrow(w7_var_df) +nrow(w8_var_df) 

#' Tidyverse method to join the dataframes separated by wave
#' #' https://www.statology.org/merge-multiple-data-frames-in-r/
#install.packages("tidyverse")
library(tidyverse)
library(data.table)
library(magrittr)
library(purrr)
library(dplyr)

#' Put all separated data frames into list
df_list <- list( var_code_df, w2_var_df, w3_var_df, w4_var_df, w5_var_df, w6_var_df, w7_var_df, w8_var_df )      
#' Merge all data frames together
merged_df = df_list %>% reduce( full_join , by = 'variables_code' )
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined_2023_03_20/")
#saveRDS( merged_df , "wave_results_analysis_inc_new_pango_stat_reshaped.rds")
saveRDS( merged_df , "wave_results_analysis_reshaped.rds")
rm(w2_var_df,w3_var_df,w4_var_df,w5_var_df,w6_var_df,w7_var_df,w8_var_df)
wave_results_analysis_reshaped = readRDS("wave_results_analysis_reshaped.rds")


#' #' Create version of wave_results_analysis_reshaped.rds with number of false positives before and after earliest true positive and also count the number of true positives
#' ** see false_positive_split.R for code **
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/false_positive_split_2023_03/" )
fp_split_df = readRDS( "fp_split_df.rds" )

wave_results_analysis_reshaped_fp = merge(  x = wave_results_analysis_reshaped 
                                          , y = fp_split_df 
                                          , by = c(  "variables_code"
                                                   , "w2_earliest_tp_EWS_date"
                                                   , "w3_earliest_tp_EWS_date"
                                                   , "w4_earliest_tp_EWS_date"
                                                   , "w5_earliest_tp_EWS_date"
                                                   , "w6_earliest_tp_EWS_date"
                                                   , "w7_earliest_tp_EWS_date"
                                                   , "w8_earliest_tp_EWS_date") 
                                          , all = TRUE
                                          )

setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined_2023_03_20/" )
saveRDS( wave_results_analysis_reshaped_fp , "wave_results_analysis_reshaped_fp.rds" ) #' same as wave_results_analysis_reshaped but with false positives split by before and after earliest true positive

#' Check that <false positive EWS before earliest true positive EWS> + <fp EWS after earliest tp EWS> = <total fp>
#' But note that the sum will be zero if there are no tp EWS
#' Filter out rows where there is no tp EWS
#' Wave 2 
tp_only = subset(wave_results_analysis_reshaped_fp , wave_results_analysis_reshaped_fp$w2_n_tp != 0 & !is.na(wave_results_analysis_reshaped_fp$w2_n_tp) )
sum( tp_only$w2_EWS_false_positive_n )
sum( as.numeric( na.omit(tp_only$w2_n_fp_before_tp ) ) ) + sum( as.numeric( na.omit(tp_only$w2_n_fp_after_tp) ) )
#' Wave 3
tp_only = subset(wave_results_analysis_reshaped_fp , wave_results_analysis_reshaped_fp$w3_n_tp != 0 & !is.na(wave_results_analysis_reshaped_fp$w3_n_tp) )
sum( tp_only$w3_EWS_false_positive_n )
sum( na.omit(tp_only$w3_n_fp_before_tp ) ) + sum( na.omit(tp_only$w3_n_fp_after_tp) )
#' Wave 4
tp_only = subset(wave_results_analysis_reshaped_fp , wave_results_analysis_reshaped_fp$w4_n_tp != 0 & !is.na(wave_results_analysis_reshaped_fp$w4_n_tp) )
sum( tp_only$w4_EWS_false_positive_n )
sum( na.omit(tp_only$w4_n_fp_before_tp ) ) + sum( na.omit(tp_only$w4_n_fp_after_tp) )
#' Wave 5
tp_only = subset(wave_results_analysis_reshaped_fp , wave_results_analysis_reshaped_fp$w5_n_tp != 0 & !is.na(wave_results_analysis_reshaped_fp$w5_n_tp) )
sum( tp_only$w5_EWS_false_positive_n )
sum( na.omit(tp_only$w5_n_fp_before_tp ) ) + sum( na.omit(tp_only$w5_n_fp_after_tp) )
#' Wave 6
tp_only = subset(wave_results_analysis_reshaped_fp , wave_results_analysis_reshaped_fp$w6_n_tp != 0 & !is.na(wave_results_analysis_reshaped_fp$w6_n_tp) )
sum( tp_only$w6_EWS_false_positive_n )
sum( na.omit(tp_only$w6_n_fp_before_tp ) ) + sum( na.omit(tp_only$w6_n_fp_after_tp) )
#' Wave 7
tp_only = subset(wave_results_analysis_reshaped_fp , wave_results_analysis_reshaped_fp$w7_n_tp != 0 & !is.na(wave_results_analysis_reshaped_fp$w7_n_tp) )
sum( tp_only$w7_EWS_false_positive_n )
sum( na.omit(tp_only$w7_n_fp_before_tp ) ) + sum( na.omit(tp_only$w7_n_fp_after_tp) )
#' Wave 8
tp_only = subset(wave_results_analysis_reshaped_fp , wave_results_analysis_reshaped_fp$w8_n_tp != 0 & !is.na(wave_results_analysis_reshaped_fp$w8_n_tp) )
sum( tp_only$w8_EWS_false_positive_n )
sum( na.omit(tp_only$w8_n_fp_before_tp ) ) + sum( na.omit(tp_only$w8_n_fp_after_tp) )
rm(tp_only)


##########
#' Interrogating merged data

#' Variables where there are no false positives at all
View( subset( merged_df 
              , merged_df[,11] == 0 & merged_df[,14] == 0 & merged_df[,17] == 0 & merged_df[,20] == 0 & merged_df[,23] == 0 & merged_df[,26] == 0 & merged_df[,29] == 0
            )
    )

#' Variables where there are no false positives for B.1.177, Alpha, Delta 1, Omicron 1 and Omicron 2 (these are the waves driven by genomic changes)
View( subset( merged_df 
              , merged_df[,11] == 0 & merged_df[,14] == 0 & merged_df[,17] == 0 & merged_df[,26] == 0 & merged_df[,29] == 0
)
)
no_fp_main_df = subset( merged_df 
                  , merged_df[,11] == 0 & merged_df[,14] == 0 & merged_df[,17] == 0 & merged_df[,26] == 0 & merged_df[,29] == 0 )
View(no_fp_main_df[,c(2,3,4,5,6,7,8,10,13,16,19,22,25,28)])
no_fp_main_2_df = no_fp_main_df[,c(2,3,4,5,6,7,8,10,13,16,19,22,25,28)]
#no_fp_main_2_df[,seq(8,14,1)] = no_fp_main_2_df[,seq(8,14,1)] #+ 7 #' +7 represents the assumed lag in reporting hospitalisation data (DO NOT NEED TO APPLY)
unique(no_fp_main_2_df[,1]) ; unique(no_fp_main_2_df[,2]) ; unique(no_fp_main_2_df[,3]); unique(no_fp_main_2_df[,4]); unique(no_fp_main_2_df[,5]); unique(no_fp_main_2_df[,6]); sort(unique(no_fp_main_2_df[,7]))
message("Lead/lag range for B.1.177 is "     ,min(no_fp_main_2_df[, 8])," to ", max(no_fp_main_2_df[, 8]) );                   
message("Lead/lag range for Alpha is "       ,min(no_fp_main_2_df[, 9])," to ", max(no_fp_main_2_df[, 9]) );                   
message("Lead/lag range for Delta 1 is "     ,min(no_fp_main_2_df[,10])," to ", max(no_fp_main_2_df[,10]) );                   
message("Lead/lag range for Delta 2 is "     ,min(na.omit(no_fp_main_2_df[,11]))," to ", max(na.omit(no_fp_main_2_df[,11])) ); 
message("Lead/lag range for Delta 3 is "     ,min(na.omit(no_fp_main_2_df[,12]))," to ", max(na.omit(no_fp_main_2_df[,12])) ); 
message("Lead/lag range for Omicron BA.1 is ",min(no_fp_main_2_df[,13])," to ", max(no_fp_main_2_df[,13]) );                   
message("Lead/lag range for Omicron BA.2 is ",min(no_fp_main_2_df[,14])," to ", max(no_fp_main_2_df[,14]) );                   

min_lead_lag = min(na.omit(no_fp_main_2_df[,c(8,9,10,11,12,13,14)]))
max_lead_lag = max(na.omit(no_fp_main_2_df[,c(8,9,10,11,12,13,14)]))
hist(no_fp_main_2_df[, 8],xlab="EWS lead(-ve) or lag(+ve)",ylab="Number of sets of parameters/cluster filters",main="B.1.177",xlim=c(-40,60),ylim=c(0,700) ) #xlim=c(min_lead_lag,max_lead_lag),ylim=c(0,700) )
hist(no_fp_main_2_df[, 9],xlab="EWS lead(-ve) or lag(+ve)",ylab="Number of sets of parameters/cluster filters",main="Alpha",xlim=c(min_lead_lag,max_lead_lag),ylim=c(0,700) )
hist(no_fp_main_2_df[,10],xlab="EWS lead(-ve) or lag(+ve)",ylab="Number of sets of parameters/cluster filters",main="Delta 1",xlim=c(min_lead_lag,max_lead_lag),ylim=c(0,700) )
hist(no_fp_main_2_df[,11],xlab="EWS lead(-ve) or lag(+ve)",ylab="Number of sets of parameters/cluster filters",main="Delta 2",xlim=c(min_lead_lag,max_lead_lag),ylim=c(0,700) )
hist(no_fp_main_2_df[,12],xlab="EWS lead(-ve) or lag(+ve)",ylab="Number of sets of parameters/cluster filters",main="Delta 3",xlim=c(min_lead_lag,max_lead_lag),ylim=c(0,700) )
hist(no_fp_main_2_df[,13],xlab="EWS lead(-ve) or lag(+ve)",ylab="Number of sets of parameters/cluster filters",main="Omicron BA.1",xlim=c(min_lead_lag,max_lead_lag),ylim=c(0,700) )
hist(no_fp_main_2_df[,14],xlab="EWS lead(-ve) or lag(+ve)",ylab="Number of sets of parameters/cluster filters",main="Omicron BA.2",xlim=c(min_lead_lag,max_lead_lag),ylim=c(0,700) )

#' Variables where there are no false positives for selected column
wave_col = 29 # 11, 14, 17,20,23,26,29
temp_df = subset( merged_df , merged_df[ , wave_col ] == 0 ) #' Filter by zero false positives for a single wave
temp_df = subset( temp_df , temp_df[ , wave_col-1 ] == min( na.omit(temp_df[ , wave_col-1 ])) ) #' Filter by the minimum lead/lag for a single wave which has already been filtered for zero false positives
View(temp_df)
sort(unique(temp_df$leading_indicator_type));unique(temp_df$TFPS_cluster_min_age);unique(temp_df$TFPS_cluster_max_age);unique(temp_df$TFPS_cluster_min_descendants);unique(temp_df$LGR_p_value_threshold);unique(temp_df$parent_sub_lgr_threshold);sort(unique(temp_df$EWS_threshold))
min(temp_df[,wave_col-1])
#' Which leading indicators not present
'%ni%' = Negate('%in%')
sort(unique(wave_results_analysis_df$leading_indicator_type))[which(sort(unique(wave_results_analysis_df$leading_indicator_type)) %ni% sort(unique(temp_df$leading_indicator_type)))]


#' Variables where there are no false positives BEFORE earliest true positive (all waves)
nrow( subset( wave_results_analysis_reshaped_fp 
              , wave_results_analysis_reshaped_fp[,2]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,2]) &
                wave_results_analysis_reshaped_fp[,3]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,3]) &
                wave_results_analysis_reshaped_fp[,4]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,4]) &
                wave_results_analysis_reshaped_fp[,5]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,5]) &
                wave_results_analysis_reshaped_fp[,6]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,6]) &
                wave_results_analysis_reshaped_fp[,7]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,7]) &
                wave_results_analysis_reshaped_fp[,8]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,8])
            )
      )
View( subset( wave_results_analysis_reshaped_fp 
              , wave_results_analysis_reshaped_fp[,2]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,2]) &
                wave_results_analysis_reshaped_fp[,3]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,3]) &
                wave_results_analysis_reshaped_fp[,4]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,4]) &
                wave_results_analysis_reshaped_fp[,5]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,5]) &
                wave_results_analysis_reshaped_fp[,6]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,6]) &
                wave_results_analysis_reshaped_fp[,7]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,7]) &
                wave_results_analysis_reshaped_fp[,8]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,8]) &
                
                wave_results_analysis_reshaped_fp[,31] == 0 &
                wave_results_analysis_reshaped_fp[,34] == 0 &
                wave_results_analysis_reshaped_fp[,37] == 0 &
                wave_results_analysis_reshaped_fp[,40] == 0 &
                wave_results_analysis_reshaped_fp[,43] == 0 &
                wave_results_analysis_reshaped_fp[,46] == 0 &
                wave_results_analysis_reshaped_fp[,49] == 0
              )
      )
#' There are 40,720 parameter sets (or variable sets) with at least one true positive EWS for all 7 waves
#' However, none of these have zero false positive EWS before the earliest true positive

#' Variables where there are no false positives BEFORE earliest true positive (all waves)
nrow( subset( wave_results_analysis_reshaped_fp 
              , #wave_results_analysis_reshaped_fp[,2]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,2]) &
                wave_results_analysis_reshaped_fp[,3]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,3]) &
                wave_results_analysis_reshaped_fp[,4]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,4]) &
                wave_results_analysis_reshaped_fp[,5]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,5]) &
                wave_results_analysis_reshaped_fp[,6]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,6]) &
                wave_results_analysis_reshaped_fp[,7]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,7]) &
                wave_results_analysis_reshaped_fp[,8]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,8])
)
)
View( subset( wave_results_analysis_reshaped_fp 
              , #wave_results_analysis_reshaped_fp[,2]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,2]) &
                wave_results_analysis_reshaped_fp[,3]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,3]) &
                wave_results_analysis_reshaped_fp[,4]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,4]) &
                wave_results_analysis_reshaped_fp[,5]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,5]) &
                wave_results_analysis_reshaped_fp[,6]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,6]) &
                wave_results_analysis_reshaped_fp[,7]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,7]) &
                wave_results_analysis_reshaped_fp[,8]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,8]) &
                
                #wave_results_analysis_reshaped_fp[,31] == 0 &
                wave_results_analysis_reshaped_fp[,34] == 0 &
                wave_results_analysis_reshaped_fp[,37] == 0 &
                wave_results_analysis_reshaped_fp[,40] == 0 &
                wave_results_analysis_reshaped_fp[,43] == 0 &
                wave_results_analysis_reshaped_fp[,46] == 0 &
                wave_results_analysis_reshaped_fp[,49] == 0
)
)
#' There are 43,059 parameter sets (or variable sets) with at least one true positive EWS for all 6/7 waves (exc. B.1.177 where there is possibly not enough data points before the wave start)
#' However, none of these have zero false positive EWS before the earliest true positive


#' Variables where there are no false positives BEFORE earliest true positive (only for B.1.177, Alpha, Delta 1, Omicron 1 and Omicron 2 (these are the waves driven by genomic changes))
nrow( subset( wave_results_analysis_reshaped_fp 
              , wave_results_analysis_reshaped_fp[,2]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,2]) &
                wave_results_analysis_reshaped_fp[,3]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,3]) &
                wave_results_analysis_reshaped_fp[,4]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,4]) &
                #wave_results_analysis_reshaped_fp[,5]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,5]) &
                #wave_results_analysis_reshaped_fp[,6]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,6]) &
                wave_results_analysis_reshaped_fp[,7]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,7]) &
                wave_results_analysis_reshaped_fp[,8]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,8])
            )
    )
View( subset( wave_results_analysis_reshaped_fp 
              , wave_results_analysis_reshaped_fp[,2]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,2]) &
                wave_results_analysis_reshaped_fp[,3]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,3]) &
                wave_results_analysis_reshaped_fp[,4]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,4]) &
                #wave_results_analysis_reshaped_fp[,5]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,5]) &
                #wave_results_analysis_reshaped_fp[,6]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,6]) &
                wave_results_analysis_reshaped_fp[,7]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,7]) &
                wave_results_analysis_reshaped_fp[,8]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,8]) &
                
                wave_results_analysis_reshaped_fp[,31] == 0 &
                wave_results_analysis_reshaped_fp[,34] == 0 &
                wave_results_analysis_reshaped_fp[,37] == 0 &
                #wave_results_analysis_reshaped_fp[,40] == 0 &
                #wave_results_analysis_reshaped_fp[,43] == 0 &
                wave_results_analysis_reshaped_fp[,46] == 0 &
                wave_results_analysis_reshaped_fp[,49] == 0
)
)
#' There are 55,828 parameter sets (or variable sets) with at least one true positive EWS for the 5 genomic variant driven waves
#' However, none of these have zero false positive EWS before the earliest true positive

#' Variables where there are no false positives BEFORE earliest true positive (only for B.1.177, Alpha, Delta 1, Omicron 1 and Omicron 2 (these are the waves driven by genomic changes))
nrow( subset( wave_results_analysis_reshaped_fp 
              , #wave_results_analysis_reshaped_fp[,2]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,2]) &
                wave_results_analysis_reshaped_fp[,3]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,3]) &
                wave_results_analysis_reshaped_fp[,4]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,4]) &
                #wave_results_analysis_reshaped_fp[,5]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,5]) &
                #wave_results_analysis_reshaped_fp[,6]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,6]) &
                wave_results_analysis_reshaped_fp[,7]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,7]) &
                wave_results_analysis_reshaped_fp[,8]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,8])
)
)
View( subset( wave_results_analysis_reshaped_fp 
              , #wave_results_analysis_reshaped_fp[,2]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,2]) &
                wave_results_analysis_reshaped_fp[,3]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,3]) &
                wave_results_analysis_reshaped_fp[,4]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,4]) &
                #wave_results_analysis_reshaped_fp[,5]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,5]) &
                #wave_results_analysis_reshaped_fp[,6]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,6]) &
                wave_results_analysis_reshaped_fp[,7]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,7]) &
                wave_results_analysis_reshaped_fp[,8]!=0 &  !is.na(wave_results_analysis_reshaped_fp[,8]) &
                
                #wave_results_analysis_reshaped_fp[,31] == 0 &
                wave_results_analysis_reshaped_fp[,34] == 0 &
                wave_results_analysis_reshaped_fp[,37] == 0 &
                #wave_results_analysis_reshaped_fp[,40] == 0 &
                #wave_results_analysis_reshaped_fp[,43] == 0 &
                wave_results_analysis_reshaped_fp[,46] == 0 &
                wave_results_analysis_reshaped_fp[,49] == 0
)
)
#' There are 59,919 parameter sets (or variable sets) with at least one true positive EWS for the 4/5 genomic variant driven waves (exc. B.1.177 where possibly not enough data points before wave start)
#' However, none of these have zero false positive EWS before the earliest true positive

#' Variables where there are no false positives BEFORE earliest true positive for selected column
wave_col = 49 # 31,34,37,40,43,46,49
temp_df = subset( wave_results_analysis_reshaped_fp , wave_results_analysis_reshaped_fp[ , wave_col ] == 0 ) #' Filter by zero false positives BEFORE earliest true positive for a single wave
temp_df = subset( temp_df , temp_df[ , (((wave_col-28)/3)*2)+14 ] == min( na.omit(temp_df[ , (((wave_col-28)/3)*2)+14 ])) ) #' Filter by the minimum lead/lag for a single wave which has already been filtered for zero false positives BEFORE earliest true positive
View(temp_df)
sort(unique(temp_df$leading_indicator_type));unique(temp_df$TFPS_cluster_min_age);unique(temp_df$TFPS_cluster_max_age);unique(temp_df$TFPS_cluster_min_descendants);unique(temp_df$LGR_p_value_threshold);unique(temp_df$parent_sub_lgr_threshold);sort(unique(temp_df$EWS_threshold))
min(temp_df[,(((wave_col-28)/3)*2)+14])
#' Which leading indicators not present
'%ni%' = Negate('%in%')
sort(unique(wave_results_analysis_df$leading_indicator_type))[which(sort(unique(wave_results_analysis_df$leading_indicator_type)) %ni% sort(unique(temp_df$leading_indicator_type)))]


###################
######' Looking at impact of parent/sub-cluster replacement on Omicron BA.1 EWS timing ######
merged_df_999 = subset( merged_df , merged_df$parent_sub_lgr_threshold == 999 )# & merged_df$w7_EWS_false_positive_n==0
merged_df_not_999_060 = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$parent_sub_lgr_threshold==0.60 ) #& merged_df$w7_EWS_false_positive_n==0
merged_df_not_999_065 = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$parent_sub_lgr_threshold==0.65 ) #& merged_df$w7_EWS_false_positive_n==0
merged_df_not_999_070 = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$parent_sub_lgr_threshold==0.70 ) #& merged_df$w7_EWS_false_positive_n==0
merged_df_not_999_075 = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$parent_sub_lgr_threshold==0.75 ) #& merged_df$w7_EWS_false_positive_n==0
merged_df_not_999_080 = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$parent_sub_lgr_threshold==0.80 ) #& merged_df$w7_EWS_false_positive_n==0
merged_df_not_999_085 = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$parent_sub_lgr_threshold==0.85 ) #& merged_df$w7_EWS_false_positive_n==0
merged_df_not_999_090 = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$parent_sub_lgr_threshold==0.90 ) #& merged_df$w7_EWS_false_positive_n==0
merged_df_not_999_095 = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$parent_sub_lgr_threshold==0.95 ) #& merged_df$w7_EWS_false_positive_n==0
merged_df_not_999_100 = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$parent_sub_lgr_threshold==1.00 ) #& merged_df$w7_EWS_false_positive_n==0

hist(  na.omit(merged_df_not_999_100$w7_lead_lag_days)
     , breaks = 20
     , col=rgb(0,0,1,0.5) # red,green,blue and transparency
     , xlim = c(0,80)
     , ylim = c(0,0.2)
     , freq = FALSE
     , xlab = "EWS lead (-ve) or lag (+ve) days"
     #, main = paste("Impact of parent/sub-cluster replacement on EWS lead/lag for Omicron BA.1 wave"))
     #, main = paste("No parent/sub-cluster replacement"))
     , main = paste("Sub-clusters replaced by parent with 100% LGR threshold"))
hist(na.omit(merged_df_not_999_060$w7_lead_lag_days), breaks = 20, col=rgb(0,1,0,0.5), add=T, freq = FALSE)
#hist(na.omit(merged_df_not_999_100$w7_lead_lag_days), breaks = 20, col=rgb(1,1,0,0.5), add=T, freq = FALSE)
legend("topright",legend=c("No parent/sub-cluster replacement","Sub-clusters replaced by parent"),col=c(rgb(0,0,1,0.5),rgb(0,1,0,0.5)),lty=1,lwd=5)
#' Compile earliest, mean and median lead/lag by parent/sub-cluster replacement LGR threshold
min(na.omit(merged_df_not_999_060$w7_lead_lag_days)) ; mean(na.omit(merged_df_not_999_060$w7_lead_lag_days)) ; median(na.omit(merged_df_not_999_060$w7_lead_lag_days))
ba1_df = data.frame()
counter = 0
for (i in c(seq(0.60,1.00,0.05),999)){
  counter = counter + 1
  temp_df = subset( merged_df , merged_df$parent_sub_lgr_threshold == i )
  ba1_df[counter , 1] = paste0(min(na.omit(temp_df$w7_lead_lag_days))," ",round(mean(na.omit(temp_df$w7_lead_lag_days)),0)," ", median(na.omit(temp_df$w7_lead_lag_days)))
  temp_df = subset( merged_df , merged_df$parent_sub_lgr_threshold == i  & merged_df$w7_EWS_false_positive_n==0 )
  ba1_df[counter , 2] =  paste0(min(na.omit(temp_df$w7_lead_lag_days))," ",round(mean(na.omit(temp_df$w7_lead_lag_days)),0)," ", median(na.omit(temp_df$w7_lead_lag_days)))
  temp_df = subset( merged_df , merged_df$parent_sub_lgr_threshold == i  & merged_df$w7_EWS_false_positive_n==1 )
  ba1_df[counter , 3] =  paste0(min(na.omit(temp_df$w7_lead_lag_days))," ",round(mean(na.omit(temp_df$w7_lead_lag_days)),0)," ", median(na.omit(temp_df$w7_lead_lag_days)))
  temp_df = subset( merged_df , merged_df$parent_sub_lgr_threshold == i  & merged_df$w7_EWS_false_positive_n==2 )
  ba1_df[counter , 4] = paste0(min(na.omit(temp_df$w7_lead_lag_days))," ",round(mean(na.omit(temp_df$w7_lead_lag_days)),0)," ", median(na.omit(temp_df$w7_lead_lag_days)))
  temp_df = subset( merged_df , merged_df$parent_sub_lgr_threshold == i  & merged_df$w7_EWS_false_positive_n==3 )
  ba1_df[counter , 5] = paste0(min(na.omit(temp_df$w7_lead_lag_days))," ",round(mean(na.omit(temp_df$w7_lead_lag_days)),0)," ", median(na.omit(temp_df$w7_lead_lag_days)))
  temp_df = subset( merged_df , merged_df$parent_sub_lgr_threshold == i  & merged_df$w7_EWS_false_positive_n==4 )
  ba1_df[counter , 6] = paste0(min(na.omit(temp_df$w7_lead_lag_days))," ",round(mean(na.omit(temp_df$w7_lead_lag_days)),0)," ", median(na.omit(temp_df$w7_lead_lag_days)))
  temp_df = subset( merged_df , merged_df$parent_sub_lgr_threshold == i  & merged_df$w7_EWS_false_positive_n==5 )
  ba1_df[counter , 7] =  paste0(min(na.omit(temp_df$w7_lead_lag_days))," ",round(mean(na.omit(temp_df$w7_lead_lag_days)),0)," ", median(na.omit(temp_df$w7_lead_lag_days)))
  temp_df = subset( merged_df , merged_df$parent_sub_lgr_threshold == i  & merged_df$w7_EWS_false_positive_n>5 & merged_df$w7_EWS_false_positive_n<=10  );
  ba1_df[counter , 8] =  paste0(min(na.omit(temp_df$w7_lead_lag_days))," ",round(mean(na.omit(temp_df$w7_lead_lag_days)),0)," ", median(na.omit(temp_df$w7_lead_lag_days)))
  temp_df = subset( merged_df , merged_df$parent_sub_lgr_threshold == i  & merged_df$w7_EWS_false_positive_n > 10 );
  ba1_df[counter , 9] =  paste0(min(na.omit(temp_df$w7_lead_lag_days))," ",round(mean(na.omit(temp_df$w7_lead_lag_days)),0)," ", median(na.omit(temp_df$w7_lead_lag_days)))
}
rownames(ba1_df) <- c("60%","65%","70%","75%","80%","85%","90%","95%","100%","None")
colnames(ba1_df) <- c("Any","0","1","2","3","4","5","6-10",">10")

#*Comparison of parent/sub-cluster replacement (no pango lineage LI) with no replacement (no pango lineage LI) and no replacement (yes pango lineage LI)**
# No parent/sub-cluster replacement and no pango lineage max(LGR) leading indicator
merged_df_999_no_pango = subset( merged_df , merged_df$parent_sub_lgr_threshold == 999 & merged_df$leading_indicator_type!="tfps_pango_lineage_max_lgr" ) # & merged_df$w7_EWS_false_positive_n==0
# Yes parent/sub-cluster replacement and no pango lineage max(LGR) leading indicator
merged_df_not_999_no_pango = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$leading_indicator_type!="tfps_pango_lineage_max_lgr" ) # & merged_df$w7_EWS_false_positive_n==0
# No parent/sub-cluster replacement and only pango lineage max(LGR) leading indicator
merged_df_999_only_pango = subset( merged_df , merged_df$parent_sub_lgr_threshold == 999 & merged_df$leading_indicator_type=="tfps_pango_lineage_max_lgr" ) # & merged_df$w7_EWS_false_positive_n==0
# Yes parent/sub-cluster replacement and only pango lineage max(LGR) leading indicator
merged_df_not_999_only_pango = subset( merged_df , merged_df$parent_sub_lgr_threshold != 999 & merged_df$leading_indicator_type=="tfps_pango_lineage_max_lgr" ) # & merged_df$w7_EWS_false_positive_n==0

hist(  na.omit(merged_df_999_no_pango$w7_lead_lag_days)
       , breaks = seq(0,80,2)#20
       , col=rgb(0,0,1,0.5) # red,green,blue and transparency
       , xlim = c(0,80)
       , ylim = c(0,0.3)
       , freq = FALSE
       , xlab = "EWS lead (-ve) or lag (+ve) days"
       #, main = paste("Impact of parent/sub-cluster replacement on EWS lead/lag for Omicron BA.1 wave"))
       #, main = paste("No parent/sub-cluster replacement"))
       #, main = paste("Sub-clusters replaced by parent with 100% LGR threshold")
       , main=paste("All TFPS leading indicators except pango lineage max(LGR)")
    )
hist(na.omit(merged_df_not_999_no_pango$w7_lead_lag_days), breaks = seq(0,80,2), col=rgb(0,1,0,0.5), add=T, freq = FALSE)
legend("topright",legend=c("No parent/sub-cluster replacement" #and no pango LI"
                           ,"Sub-clusters replaced by parent"# and no pango LI"
                           #,"No parent/sub-cluster replacement and only pango LI"
                           #,"Sub-clusters replaced by parent and only pango LI"
                           )
       ,col=c(rgb(0,0,1,0.5),rgb(0,1,0,0.5),rgb(1,1,0,0.5),rgb(1,1,1,0.5)),lty=1,lwd=5)

hist(  na.omit(merged_df_999_only_pango$w7_lead_lag_days)
       , breaks = seq(0,80,2)#20
       , col=rgb(1,0,1,0.5) # red,green,blue and transparency
       , xlim = c(0,80)
       , ylim = c(0,0.3)
       , freq = FALSE
       , xlab = "EWS lead (-ve) or lag (+ve) days"
       #, main = paste("Impact of parent/sub-cluster replacement on EWS lead/lag for Omicron BA.1 wave"))
       #, main = paste("No parent/sub-cluster replacement"))
       #, main = paste("Sub-clusters replaced by parent with 100% LGR threshold")
       , main="Only pango lineage max(LGR) leading indicator")
#hist(na.omit(merged_df_999_only_pango$w7_lead_lag_days), breaks = seq(0,80,2), col=rgb(1,1,0,0.5), freq=FALSE)#add=T, freq = FALSE)
hist(na.omit(merged_df_not_999_only_pango$w7_lead_lag_days), breaks = seq(0,80,2), col=rgb(0.5,0.5,0.5,0.5), add=T, freq = FALSE)

legend("topright",legend=c(#"No parent/sub-cluster replacement and no pango LI"
                           #,"Sub-clusters replaced by parent and no pango LI"
                           "No parent/sub-cluster replacement"# and only pango LI"
                           ,"Sub-clusters replaced by parent"# and only pango LI")
                            )
       ,col=c(rgb(1,0,1,0.5),rgb(0.5,0.5,0.5,0.5)),lty=1,lwd=5)
################################################################################

## Heat map #######################
#' Create matrix for heat map (note that the heatmap is plotted in the opposite vertical direction to how the matrix is views)
#' Number of false positives in rows
rows_n = max( na.omit( wave_results_analysis_df$ews_false_positive_n ) ) + 1 #' +1 as need to include 0 false positives
#' Earliest lead/lag in columns
cols_n = (max( na.omit( wave_results_analysis_df$lead_lag_days ) )) - (min( na.omit( wave_results_analysis_df$lead_lag_days ) )) + 1 
wave_2_matrix = matrix( data = NA , nrow = rows_n , ncol = cols_n )
wave_3_matrix = matrix( data = NA , nrow = rows_n , ncol = cols_n )
wave_4_matrix = matrix( data = NA , nrow = rows_n , ncol = cols_n )
wave_5_matrix = matrix( data = NA , nrow = rows_n , ncol = cols_n )
wave_6_matrix = matrix( data = NA , nrow = rows_n , ncol = cols_n )
wave_7_matrix = matrix( data = NA , nrow = rows_n , ncol = cols_n )
wave_8_matrix = matrix( data = NA , nrow = rows_n , ncol = cols_n )

for ( r in 1:rows_n ){
  for ( c in 1:cols_n ){
    #' adjusted search values to place in row r (EWS lead/lag days) and col c (EWS false positive n)
    r_adj =  r - 1 #' -1 adjusts for zero false positive EWS being in the 1st row
    c_adj = c + ( min( na.omit( wave_results_analysis_df$lead_lag_days ) ) ) -1 #' adjusts for range of lead/lag days not beginning at 1 (1st column in dataframe)
    wave_2_matrix[r,c] = nrow( subset( merged_df , merged_df[,10] == c_adj & merged_df[,11] == r_adj ) )
    wave_3_matrix[r,c] = nrow( subset( merged_df , merged_df[,13] == c_adj & merged_df[,14] == r_adj ) )
    wave_4_matrix[r,c] = nrow( subset( merged_df , merged_df[,16] == c_adj & merged_df[,17] == r_adj ) )
    wave_5_matrix[r,c] = nrow( subset( merged_df , merged_df[,19] == c_adj & merged_df[,20] == r_adj ) )
    wave_6_matrix[r,c] = nrow( subset( merged_df , merged_df[,22] == c_adj & merged_df[,23] == r_adj ) )
    wave_7_matrix[r,c] = nrow( subset( merged_df , merged_df[,25] == c_adj & merged_df[,26] == r_adj ) )
    wave_8_matrix[r,c] = nrow( subset( merged_df , merged_df[,28] == c_adj & merged_df[,29] == r_adj ) )
    message("EWS lead/lag: ",r_adj," EWS false positive: ",c_adj)
  }
}

#' Above for loop re-written for speed
#'  Using data table is faster than dataframe. Also pare down the dataframe for use in loop
merged_dt = as.data.table(merged_df)
w2_dt = merged_dt[,10:11]; w3_dt = merged_dt[,13:14]; w4_dt = merged_dt[,16:17] ; w5_dt = merged_dt[,19:20]
w6_dt = merged_dt[,22:23]; w7_dt = merged_dt[,25:26]; w8_dt = merged_dt[,28:29]
rm(merged_dt)
for ( r in 1:rows_n ){
  for ( c in 1:cols_n ){
    #' adjusted search values to place in row r (EWS lead/lag days) and col c (EWS false positive n)
    r_adj =  r - 1 #' -1 adjusts for zero false positive EWS being in the 1st row
    c_adj = c + ( min( na.omit( wave_results_analysis_df$lead_lag_days ) ) ) -1 #' adjusts for range of lead/lag days not beginning at 1 (1st column in dataframe)
    wave_2_matrix[r,c] = nrow( w2_dt[ w2_lead_lag_days+7 == c_adj & w2_EWS_false_positive_n == r_adj ] )
    wave_3_matrix[r,c] = nrow( w3_dt[ w3_lead_lag_days+7 == c_adj & w3_EWS_false_positive_n == r_adj ] )
    wave_4_matrix[r,c] = nrow( w4_dt[ w4_lead_lag_days+7 == c_adj & w4_EWS_false_positive_n == r_adj ] )
    wave_5_matrix[r,c] = nrow( w5_dt[ w5_lead_lag_days+7 == c_adj & w5_EWS_false_positive_n == r_adj ] )
    wave_6_matrix[r,c] = nrow( w6_dt[ w6_lead_lag_days+7 == c_adj & w6_EWS_false_positive_n == r_adj ] )
    wave_7_matrix[r,c] = nrow( w7_dt[ w7_lead_lag_days+7 == c_adj & w7_EWS_false_positive_n == r_adj ] )
    wave_8_matrix[r,c] = nrow( w8_dt[ w8_lead_lag_days+7 == c_adj & w8_EWS_false_positive_n == r_adj ] )
    message("EWS lead/lag: ",c_adj," EWS false positive: ",r_adj)
  }
}
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/wave_matrices_for_heatmap/")
saveRDS(wave_2_matrix , "wave_2_matrix.rds") ; saveRDS(wave_3_matrix , "wave_3_matrix.rds") ; saveRDS(wave_4_matrix , "wave_4_matrix.rds") ; saveRDS(wave_5_matrix , "wave_5_matrix.rds") ; saveRDS(wave_6_matrix , "wave_6_matrix.rds") ; saveRDS(wave_7_matrix , "wave_7_matrix.rds") ; saveRDS(wave_8_matrix , "wave_8_matrix.rds")

wave_2_matrix_perc = ( wave_2_matrix / sum(wave_2_matrix) ) * 100
wave_3_matrix_perc = ( wave_3_matrix / sum(wave_3_matrix) ) * 100
wave_4_matrix_perc = ( wave_4_matrix / sum(wave_4_matrix) ) * 100
wave_5_matrix_perc = ( wave_5_matrix / sum(wave_5_matrix) ) * 100
wave_6_matrix_perc = ( wave_6_matrix / sum(wave_6_matrix) ) * 100
wave_7_matrix_perc = ( wave_7_matrix / sum(wave_7_matrix) ) * 100
wave_8_matrix_perc = ( wave_8_matrix / sum(wave_8_matrix) ) * 100

#'**WORK IN PROGRESS**
#' Possibly faster to use an array
ar = array( seq(1,7*89*101,1) , dim = c( cols_n , row_n , 7 ) ) #' dimensions number of lead/lag days, number of EWS false positives, SC2 UK infection wave e.g. 1 = B.1.177
#' and a list of combinations with lapply
data.table::CJ( w2_dt[ w2_lead_lag_days] , w2_dt[ w2_EWS_false_positive_n ] )


#' test heatmap
test_matrix = matrix(data=NA,nrow=2,ncol=3)
test_matrix[1,1]=1
test_matrix[1,2]=2
test_matrix[1,3]=0
test_matrix[2,1]=4
test_matrix[2,2]=5
test_matrix[2,3]=2000
heatmap(test_matrix, Rowv = NA, Colv = NA) #' Colour ranges for heatmap mean that if most values are similar and one is a lot higher then the colour of the higher value will be very similar to the next highest, even if it is much smaller
View(test_matrix)
plot_ly(z=test_matrix, type = "heatmap",x=c(1,2,3),y=c(1,2)) #' Colour ranges for heatmap show clear difference between large and small values but if there is a much larger value then the smaller values that are different cannot be distinguished

plot_ly(z=wave_2_matrix, type = "heatmap")#, x = lead_lag_labels, y = ews_false_positive_labels )
#' Change row and column labels
rownames(wave_2_matrix) = seq(0,rows_n-1,1) ; colnames(wave_2_matrix) = seq( min(na.omit(wave_results_analysis_df$lead_lag_days)) , max(na.omit(wave_results_analysis_df$lead_lag_days)) , 1 )
rownames(wave_3_matrix) = seq(0,rows_n-1,1) ; colnames(wave_3_matrix) = seq( min(na.omit(wave_results_analysis_df$lead_lag_days)) , max(na.omit(wave_results_analysis_df$lead_lag_days)) , 1 )
rownames(wave_4_matrix) = seq(0,rows_n-1,1) ; colnames(wave_4_matrix) = seq( min(na.omit(wave_results_analysis_df$lead_lag_days)) , max(na.omit(wave_results_analysis_df$lead_lag_days)) , 1 )
rownames(wave_5_matrix) = seq(0,rows_n-1,1) ; colnames(wave_5_matrix) = seq( min(na.omit(wave_results_analysis_df$lead_lag_days)) , max(na.omit(wave_results_analysis_df$lead_lag_days)) , 1 )
rownames(wave_6_matrix) = seq(0,rows_n-1,1) ; colnames(wave_6_matrix) = seq( min(na.omit(wave_results_analysis_df$lead_lag_days)) , max(na.omit(wave_results_analysis_df$lead_lag_days)) , 1 )
rownames(wave_7_matrix) = seq(0,rows_n-1,1) ; colnames(wave_7_matrix) = seq( min(na.omit(wave_results_analysis_df$lead_lag_days)) , max(na.omit(wave_results_analysis_df$lead_lag_days)) , 1 )
rownames(wave_8_matrix) = seq(0,rows_n-1,1) ; colnames(wave_8_matrix) = seq( min(na.omit(wave_results_analysis_df$lead_lag_days)) , max(na.omit(wave_results_analysis_df$lead_lag_days)) , 1 )

heatmap( wave_2_matrix_perc, Rowv = NA, Colv = NA, scale="none")#, labcol = ews_false_positive_labels , rowcol = lead_lag_labels)# ,scale = c("row"))#, x = lead_lag_labels, y = ews_false_positive_labels )
heatmap( wave_3_matrix_perc, Rowv = NA, Colv = NA)#, x = lead_lag_labels, y = ews_false_positive_labels )
heatmap( wave_4_matrix_perc, Rowv = NA, Colv = NA)#, x = lead_lag_labels, y = ews_false_positive_labels )
heatmap( wave_5_matrix_perc, Rowv = NA, Colv = NA)#, x = lead_lag_labels, y = ews_false_positive_labels )
heatmap( wave_6_matrix_perc, Rowv = NA, Colv = NA)#, x = lead_lag_labels, y = ews_false_positive_labels )
heatmap( wave_7_matrix_perc, Rowv = NA, Colv = NA)#, x = lead_lag_labels, y = ews_false_positive_labels )
heatmap( wave_8_matrix_perc, Rowv = NA, Colv = NA)#, x = lead_lag_labels, y = ews_false_positive_labels )

library(plotly)
library(RColorBrewer)
# your palette definition
heatmap_palette <- colorRampPalette(c("darkblue", "blue", "lightblue1","green","yellow", "red", "darkred"))
heatmap_palette <- colorRampPalette(c("white","yellow", "orange", "red", "darkred"))

lead_lag_labels = seq( min(na.omit(wave_results_analysis_df$lead_lag_days))+7 , max(na.omit(wave_results_analysis_df$lead_lag_days))+7 , 1 )
ews_false_positive_labels = seq(0,rows_n-1,1)
plot_ly(z = wave_2_matrix , colors = heatmap_palette(100000) , type = "heatmap", x = lead_lag_labels, y = ews_false_positive_labels ); 
plot_ly(z = wave_2_matrix , type = "heatmap", x = lead_lag_labels, y = ews_false_positive_labels )
plot_ly(z = wave_2_matrix_perc , type = "heatmap", x = lead_lag_labels, y = ews_false_positive_labels )
plot_ly(z = wave_2_matrix_perc[1:10,] , type = "heatmap", x = lead_lag_labels, y = ews_false_positive_labels )
plot_ly(z = wave_3_matrix_perc , type = "heatmap", x = lead_lag_labels, y = ews_false_positive_labels )
plot_ly(z = wave_4_matrix_perc , type = "heatmap", x = lead_lag_labels, y = ews_false_positive_labels )
plot_ly(z = wave_5_matrix_perc , type = "heatmap", x = lead_lag_labels, y = ews_false_positive_labels )
plot_ly(z = wave_6_matrix_perc , type = "heatmap", x = lead_lag_labels, y = ews_false_positive_labels )
plot_ly(z = wave_7_matrix_perc , type = "heatmap", x = lead_lag_labels, y = ews_false_positive_labels )
plot_ly(z = wave_8_matrix_perc , type = "heatmap", x = lead_lag_labels, y = ews_false_positive_labels )
  
#wave_2_table = table( merged_df[,11] , merged_df[,10] +7 )
#wave_3_table = table( merged_df[,14] , merged_df[,13] +7 )
#wave_4_table = table( merged_df[,17] , merged_df[,16] +7 )
#wave_5_table = table( unique(merged_df[,20]) , merged_df[,19] +7 )
#wave_6_table = table( merged_df[,23] , merged_df[,22] +7 )
#wave_7_table = table( merged_df[,26] , merged_df[,25] +7 )
#wave_8_table = table( merged_df[,29] , merged_df[,28] +7 )

#heatmap( wave_2_table, Rowv = NA, Colv = NA )
#heatmap( wave_3_table, Rowv = NA, Colv = NA )
#heatmap( wave_4_table, Rowv = NA, Colv = NA )
#heatmap( wave_5_table, Rowv = NA, Colv = NA )
#heatmap( wave_6_table, Rowv = NA, Colv = NA )
#heatmap( wave_7_table, Rowv = NA, Colv = NA )
#heatmap( wave_8_table, Rowv = NA, Colv = NA )

#View(wave_2_table)
for (i in 1:nrow(wave_2_table) ){
  wave_2_matrix[ wave_2_table[i,1] , wave_2_table[i,2] ] = wave_2_table[i,3]
}

################################################################################

#' Plot point graphs
#' Wave 2 - B.1.177
plot( merged_df[,10] , merged_df[,11]
      , xlim = c( min( na.omit( wave_results_analysis_df$lead_lag_days ) ) , max( na.omit( wave_results_analysis_df$lead_lag_days ) ) )
      , ylim = c(0, max( na.omit( wave_results_analysis_df$ews_false_positive_n ) ) )
      , xlab = "EWS lead(-ve) or lag(+ve) (days)"
      , ylab = "Number of false positive EWS"
      #, main = "B.1.177"
)
abline(v=0)

#' Wave 3 - Alpha
plot( merged_df[,13] , merged_df[,14]
      , xlim = c( min( na.omit( wave_results_analysis_df$lead_lag_days ) ) , max( na.omit( wave_results_analysis_df$lead_lag_days ) ) )
      , ylim = c(0, max( na.omit( wave_results_analysis_df$ews_false_positive_n ) ) )
      , xlab = "EWS lead(-ve) or lag(+ve) (days)"
      , ylab = "Number of false positive EWS"
      #, main = "Alpha"
)
abline(v=0)
#' Wave 4 - Delta 1
plot( merged_df[,16] , merged_df[,17]
      , xlim = c( min( na.omit( wave_results_analysis_df$lead_lag_days ) ) , max( na.omit( wave_results_analysis_df$lead_lag_days ) ) )
      , ylim = c(0, max( na.omit( wave_results_analysis_df$ews_false_positive_n ) ) )
      , xlab = "EWS lead(-ve) or lag(+ve) (days)"
      , ylab = "Number of false positive EWS"
      #, main = "Delta 1"
)
abline(v=0)
#' Wave 5 - Delta 2
plot( merged_df[,19] , merged_df[,20]
      , xlim = c( min( na.omit( wave_results_analysis_df$lead_lag_days ) ) , max( na.omit( wave_results_analysis_df$lead_lag_days ) ) )
      , ylim = c(0, max( na.omit( wave_results_analysis_df$ews_false_positive_n ) ) )
      , xlab = "EWS lead(-ve) or lag(+ve) (days)"
      , ylab = "Number of false positive EWS"
      #, main = "Delta 2"
)
abline(v=0)
#' Wave 6 - Delta 3
plot( merged_df[,22] , merged_df[,23]
      , xlim = c( min( na.omit( wave_results_analysis_df$lead_lag_days ) ) , max( na.omit( wave_results_analysis_df$lead_lag_days ) ) )
      , ylim = c(0, max( na.omit( wave_results_analysis_df$ews_false_positive_n ) ) )
      , xlab = "EWS lead(-ve) or lag(+ve) (days)"
      , ylab = "Number of false positive EWS"
      #, main = "Delta 3"
)
abline(v=0)
#' Wave 7 - Omicron BA.1
plot( merged_df[,25] , merged_df[,26]
      , xlim = c( min( na.omit( wave_results_analysis_df$lead_lag_days ) ) , max( na.omit( wave_results_analysis_df$lead_lag_days ) ) )
      , ylim = c(0, max( na.omit( wave_results_analysis_df$ews_false_positive_n ) ) )
      , xlab = "EWS lead(-ve) or lag(+ve) (days)"
      , ylab = "Number of false positive EWS"
      #, main = "Omicron BA.1"
)
abline(v=0)
#' Wave 8 - Omicron BA.2
plot( merged_df[,28] , merged_df[,29]
      , xlim = c( min( na.omit( wave_results_analysis_df$lead_lag_days ) ) , max( na.omit( wave_results_analysis_df$lead_lag_days ) ) )
      , ylim = c(0, max( na.omit( wave_results_analysis_df$ews_false_positive_n ) ) )
      , xlab = "EWS lead(-ve) or lag(+ve) (days)"
      , ylab = "Number of false positive EWS"
      #, main = "Omicron BA.2"
)
abline(v=0)

#############

#' cycle through all variables and plot fit of EWS lead or lag against number of false positive EWS
plot( 0,0 
      , xlim = c( min( na.omit( wave_results_analysis_df$lead_lag_days ) ) , max( na.omit( wave_results_analysis_df$lead_lag_days ) ) )
      , ylim = c(0, max( na.omit( wave_results_analysis_df$ews_false_positive_n ) ) )
      , xlab = "EWS lead(-ve) or lag(+ve) (days)"
      , ylab = "Number of false positive EWS"
)
counter = 0
for ( wave_no in na.omit( unique( wave_results_analysis_df$wave_n ) ) ){
  for ( lead_ind in na.omit( unique( wave_results_analysis_df$leading_indicator_type ) ) ){
    for ( min_age in na.omit( unique( wave_results_analysis_df$TFPS_cluster_min_age ) ) ){
      for ( max_age in na.omit( unique( wave_results_analysis_df$TFPS_cluster_max_age ) ) ){
        for ( min_desc in na.omit( unique( wave_results_analysis_df$TFPS_cluster_min_descendants ) ) ){
          for ( p_val_th in na.omit( unique( wave_results_analysis_df$LGR_p_value_threshold ) ) ){
            for ( lgr_th in na.omit( unique( wave_results_analysis_df$parent_sub_lgr_threshold ) ) ){
              #' Loop through EWS threshold levels, which represent the classifier boundaries
              #for ( ews_th in unique( wave_results_analysis_df$EWS_threshold ) ){
                message( wave_no," ",lead_ind," ", min_age," ",max_age," ",min_desc," ",p_val_th, " ", lgr_th)#," ",ews_th )
                #' Filter by leading indicator, scan variables, cluster filter (LGR p-value) and EWS threshold
                #' Plot relationship between number of false positives and EWS lead/lag
                df = subset( wave_results_analysis_df , wave_results_analysis_df$wave_n == wave_no &
                             wave_results_analysis_df$leading_indicator_type == lead_ind &
                             wave_results_analysis_df$TFPS_cluster_min_age == min_age &
                             wave_results_analysis_df$TFPS_cluster_max_age == max_age &
                             wave_results_analysis_df$TFPS_cluster_min_descendants == min_desc &
                             wave_results_analysis_df$LGR_p_value_threshold == p_val_th &
                             wave_results_analysis_df$parent_sub_lgr_threshold == lgr_th
                            )
                if ( sum( !is.na( df$lead_lag_days ) ) > 0 ){
                  df = data.table::setorder( df , "lead_lag_days" )
                  #counter=counter+1
                  x = df$lead_lag_days
                  y = df$ews_false_positive_n
                  #lines( x , y , col = counter , pch = 16 )
                  fit1 <- lm( y ~ x , data=df )
                  x_axis = seq(-100,100,1)
                  lines( x_axis , predict( fit1 , data.frame( x = x_axis ) ) , col = wave_no )#counter )
                }
              #}
            }
          }
        }
      }
    }
  }
}

#' Select particular variables and cycle through the EWS thresholds

wave_no = 4
lead_ind = "tfps_vlgr_samp" #substring( filenames_prefix[ 1 ], 1 , nchar( filenames_prefix[ 1 ] ) - 1 ) #"tfps_vlgr_samp"
min_age = "14"
max_age = "84"
min_desc = "100"
p_val_th = 10000#0.01 #0.05 #
ps_lgr_th = 1 #"false" #0.6
false_pos = 0

df = subset( wave_results_analysis_df 
             , wave_results_analysis_df$wave_n                       == wave_no &
               wave_results_analysis_df$leading_indicator_type       == lead_ind &
               wave_results_analysis_df$TFPS_cluster_min_age         == min_age &
               wave_results_analysis_df$TFPS_cluster_max_age         == max_age &
               wave_results_analysis_df$TFPS_cluster_min_descendants == min_desc &
               wave_results_analysis_df$LGR_p_value_threshold        == p_val_th &
               wave_results_analysis_df$parent_sub_LGR_threshold     == ps_lgr_th 
               #wave_results_analysis_df$ews_false_positive_n == false_pos 
             #wave_results_analysis_df$lead_lag_days == 
)
df = data.table::setorder( df , "lead_lag_days" )
View(df)

x = df$lead_lag_days
y = df$ews_false_positive_n
#' Filter to get the minimum number of false positives for each number of lead / lag days
#x = unique( df$lead_lag_days )
#y = lapply( subset( df , df$lead_lag_days == unique( df$lead_lag_days ) ) , min )

plot( x
      , y 
      , xlab = "Earliest EWS lead(-ve) or lag(+ve) (days)"
      , ylab = "Number of false positive EWS"
      , typ="o"
      , main = lead_ind
)
fit1 <- lm( y ~ x , data=df )
x_axis = seq(-10,40,1)
lines(x_axis, predict(fit1, data.frame(x=x_axis)), col='green')

plot( df$EWS_threshold , df$lead_lag_days , typ = "l" , ylim = c(0,50), ylab = "", xlab="EWS threshold level ('robust' Z-score)")
lines( df$EWS_threshold , df$ews_false_positive_n , col = "red")
legend("topleft",legend=c("Earliest true EWS lead (-ve) or lag (+ve)","Number of false positive EWS"),col=c("black","red"),lty=1,lwd=1)


#' Plot of lead/lag days and wave number
plot(df$lead_lag_days,df$wave_n,xlab="lead or lag days",ylab="wave number",main=lead_ind)




#############
#' Rating the ability of different parameter sets in generating EWS
#' using total number of false positives
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined_2023_03_20/')
#wave_results_analysis_reshaped_df = readRDS("wave_results_analysis_reshaped.rds") #' note the addition of _df once file read
wave_results_analysis_reshaped_fp_df = readRDS("wave_results_analysis_reshaped_fp.rds") #' note the addition of _df once file read
rm( wave_results_analysis_reshaped , wave_results_analysis_reshaped_df )

#' Replace NA in wx_n_fp_before_tp columns with 0 (have checked that all NAs relate to parameter sets where the total number of false positives is zero)
wave_results_analysis_reshaped_fp_df$w2_n_fp_before_tp[ is.na( wave_results_analysis_reshaped_fp_df$w2_n_fp_before_tp ) ] <- 0
wave_results_analysis_reshaped_fp_df$w3_n_fp_before_tp[ is.na( wave_results_analysis_reshaped_fp_df$w3_n_fp_before_tp ) ] <- 0
wave_results_analysis_reshaped_fp_df$w4_n_fp_before_tp[ is.na( wave_results_analysis_reshaped_fp_df$w4_n_fp_before_tp ) ] <- 0
wave_results_analysis_reshaped_fp_df$w5_n_fp_before_tp[ is.na( wave_results_analysis_reshaped_fp_df$w5_n_fp_before_tp ) ] <- 0
wave_results_analysis_reshaped_fp_df$w6_n_fp_before_tp[ is.na( wave_results_analysis_reshaped_fp_df$w6_n_fp_before_tp ) ] <- 0
wave_results_analysis_reshaped_fp_df$w7_n_fp_before_tp[ is.na( wave_results_analysis_reshaped_fp_df$w7_n_fp_before_tp ) ] <- 0
wave_results_analysis_reshaped_fp_df$w8_n_fp_before_tp[ is.na( wave_results_analysis_reshaped_fp_df$w8_n_fp_before_tp ) ] <- 0

#' Filter for variable/parameter sets that produce a true positive for all waves
wave_results_analysis_reshaped_fp_df = subset( wave_results_analysis_reshaped_fp_df
                                            , !is.na( wave_results_analysis_reshaped_fp_df$w2_earliest_tp_EWS_date ) &
                                              !is.na( wave_results_analysis_reshaped_fp_df$w3_earliest_tp_EWS_date ) &
                                              !is.na( wave_results_analysis_reshaped_fp_df$w4_earliest_tp_EWS_date ) &
                                              !is.na( wave_results_analysis_reshaped_fp_df$w5_earliest_tp_EWS_date ) &
                                              !is.na( wave_results_analysis_reshaped_fp_df$w6_earliest_tp_EWS_date ) &
                                              !is.na( wave_results_analysis_reshaped_fp_df$w7_earliest_tp_EWS_date ) &
                                              !is.na( wave_results_analysis_reshaped_fp_df$w8_earliest_tp_EWS_date )
)
#' 40,720 parameter sets produce a true positive for all waves


#' Create sums of lead times and false positives (total and only before earliest true positive) for all waves and 
#' also 'key' waves (those driven by new genomic variants - Alpha, Delta 1, BA.1 and BA.2)
#' Can then rank by the different totals
lead_time_total_all = wave_results_analysis_reshaped_fp_df$w2_lead_lag_days +
                      wave_results_analysis_reshaped_fp_df$w3_lead_lag_days +
                      wave_results_analysis_reshaped_fp_df$w4_lead_lag_days +
                      wave_results_analysis_reshaped_fp_df$w5_lead_lag_days +
                      wave_results_analysis_reshaped_fp_df$w6_lead_lag_days +
                      wave_results_analysis_reshaped_fp_df$w7_lead_lag_days +
                      wave_results_analysis_reshaped_fp_df$w8_lead_lag_days
lead_time_total_key = #wave_results_analysis_reshaped_fp_df$w2_lead_lag_days +
                      wave_results_analysis_reshaped_fp_df$w3_lead_lag_days +
                      wave_results_analysis_reshaped_fp_df$w4_lead_lag_days +
                      #wave_results_analysis_reshaped_fp_df$w5_lead_lag_days +
                      #wave_results_analysis_reshaped_fp_df$w6_lead_lag_days +
                      wave_results_analysis_reshaped_fp_df$w7_lead_lag_days +
                      wave_results_analysis_reshaped_fp_df$w8_lead_lag_days

false_positive_total_all =  wave_results_analysis_reshaped_fp_df$w2_EWS_false_positive_n +
                            wave_results_analysis_reshaped_fp_df$w3_EWS_false_positive_n +
                            wave_results_analysis_reshaped_fp_df$w4_EWS_false_positive_n +
                            wave_results_analysis_reshaped_fp_df$w5_EWS_false_positive_n +
                            wave_results_analysis_reshaped_fp_df$w6_EWS_false_positive_n +
                            wave_results_analysis_reshaped_fp_df$w7_EWS_false_positive_n +
                            wave_results_analysis_reshaped_fp_df$w8_EWS_false_positive_n
false_positive_total_key =  #wave_results_analysis_reshaped_fp_df$w2_EWS_false_positive_n +
                            wave_results_analysis_reshaped_fp_df$w3_EWS_false_positive_n +
                            wave_results_analysis_reshaped_fp_df$w4_EWS_false_positive_n +
                            #wave_results_analysis_reshaped_fp_df$w5_EWS_false_positive_n +
                            #wave_results_analysis_reshaped_fp_df$w6_EWS_false_positive_n +
                            wave_results_analysis_reshaped_fp_df$w7_EWS_false_positive_n +
                            wave_results_analysis_reshaped_fp_df$w8_EWS_false_positive_n
#' b4tp = number of false positives before the earliest true positive in each wave
false_positive_b4tp_total_all = wave_results_analysis_reshaped_fp_df$w2_n_fp_before_tp +
                                wave_results_analysis_reshaped_fp_df$w3_n_fp_before_tp +
                                wave_results_analysis_reshaped_fp_df$w4_n_fp_before_tp +
                                wave_results_analysis_reshaped_fp_df$w5_n_fp_before_tp +
                                wave_results_analysis_reshaped_fp_df$w6_n_fp_before_tp +
                                wave_results_analysis_reshaped_fp_df$w7_n_fp_before_tp +
                                wave_results_analysis_reshaped_fp_df$w8_n_fp_before_tp
false_positive_b4tp_total_key = #wave_results_analysis_reshaped_fp_df$w2_n_fp_before_tp +
                                wave_results_analysis_reshaped_fp_df$w3_n_fp_before_tp +
                                wave_results_analysis_reshaped_fp_df$w4_n_fp_before_tp +
                                #wave_results_analysis_reshaped_fp_df$w5_n_fp_before_tp +
                                #wave_results_analysis_reshaped_fp_df$w6_n_fp_before_tp +
                                wave_results_analysis_reshaped_fp_df$w7_n_fp_before_tp +
                                wave_results_analysis_reshaped_fp_df$w8_n_fp_before_tp

total_score_all        = lead_time_total_all + false_positive_total_all
total_score_key        = lead_time_total_key + false_positive_total_key
total_score_fpb4tp_all = lead_time_total_all + false_positive_b4tp_total_all
total_score_fpb4tp_key = lead_time_total_key + false_positive_b4tp_total_key

wave_results_analysis_reshaped_fp_df = cbind( wave_results_analysis_reshaped_fp_df 
                                              , lead_time_total_all , lead_time_total_key
                                              , false_positive_total_all , false_positive_total_key, false_positive_b4tp_total_all, false_positive_b4tp_total_key
                                              , total_score_all , total_score_key, total_score_fpb4tp_all, total_score_fpb4tp_key )
rm( lead_time_total_all , lead_time_total_key , false_positive_total_all , false_positive_total_key, false_positive_b4tp_total_all, false_positive_b4tp_total_key    , total_score_all , total_score_key, total_score_fpb4tp_all, total_score_fpb4tp_key )

#' Add indicators for maximum number of false positive EWS allowed per wave so can filter. By All & key waves and total FP and FP before TP.
fp_limits = c(0,2,5,10)
for (m in fp_limits){
  temp = ifelse(  wave_results_analysis_reshaped_fp_df$w2_EWS_false_positive_n <= m &
                  wave_results_analysis_reshaped_fp_df$w3_EWS_false_positive_n <= m &
                  wave_results_analysis_reshaped_fp_df$w4_EWS_false_positive_n <= m &
                  wave_results_analysis_reshaped_fp_df$w5_EWS_false_positive_n <= m &
                  wave_results_analysis_reshaped_fp_df$w6_EWS_false_positive_n <= m &
                  wave_results_analysis_reshaped_fp_df$w7_EWS_false_positive_n <= m &
                  wave_results_analysis_reshaped_fp_df$w8_EWS_false_positive_n <= m
                  , "y","n")
  new_name = paste0("fp_limit_",m,"_all")
  assign( new_name , temp )
}
table(fp_limit_0_all);table(fp_limit_2_all);table(fp_limit_5_all);table(fp_limit_10_all)
for (m in fp_limits){
  temp = ifelse(    wave_results_analysis_reshaped_fp_df$w3_EWS_false_positive_n <= m &
                    wave_results_analysis_reshaped_fp_df$w4_EWS_false_positive_n <= m &
                    wave_results_analysis_reshaped_fp_df$w7_EWS_false_positive_n <= m &
                    wave_results_analysis_reshaped_fp_df$w8_EWS_false_positive_n <= m
                  , "y","n")
  new_name = paste0("fp_limit_",m,"_key")
  assign( new_name , temp )
}
table(fp_limit_0_key);table(fp_limit_2_key);table(fp_limit_5_key);table(fp_limit_10_key)
for (m in fp_limits){
  temp = ifelse(  wave_results_analysis_reshaped_fp_df$w2_n_fp_before_tp <= m &
                    wave_results_analysis_reshaped_fp_df$w3_n_fp_before_tp <= m &
                    wave_results_analysis_reshaped_fp_df$w4_n_fp_before_tp <= m &
                    wave_results_analysis_reshaped_fp_df$w5_n_fp_before_tp <= m &
                    wave_results_analysis_reshaped_fp_df$w6_n_fp_before_tp <= m &
                    wave_results_analysis_reshaped_fp_df$w7_n_fp_before_tp <= m &
                    wave_results_analysis_reshaped_fp_df$w8_n_fp_before_tp <= m
                  , "y","n")
  new_name = paste0("fp_b4tp_limit_",m,"_all")
  assign( new_name , temp )
}
table(fp_b4tp_limit_0_all);table(fp_b4tp_limit_2_all);table(fp_b4tp_limit_5_all);table(fp_b4tp_limit_10_all)
for (m in fp_limits){
  temp = ifelse(    wave_results_analysis_reshaped_fp_df$w3_n_fp_before_tp <= m &
                    wave_results_analysis_reshaped_fp_df$w4_n_fp_before_tp <= m &
                    wave_results_analysis_reshaped_fp_df$w7_n_fp_before_tp <= m &
                    wave_results_analysis_reshaped_fp_df$w8_n_fp_before_tp <= m
                  , "y","n")
  new_name = paste0("fp_b4tp_limit_",m,"_key")
  assign( new_name , temp )
}
table(fp_b4tp_limit_0_key);table(fp_b4tp_limit_2_key);table(fp_b4tp_limit_5_key);table(fp_b4tp_limit_10_key)

wave_results_analysis_reshaped_fp_df = cbind( wave_results_analysis_reshaped_fp_df 
                                             , fp_limit_0_all, fp_limit_2_all, fp_limit_5_all, fp_limit_10_all
                                             , fp_limit_0_key, fp_limit_2_key, fp_limit_5_key, fp_limit_10_key
                                             , fp_b4tp_limit_0_all, fp_b4tp_limit_2_all, fp_b4tp_limit_5_all, fp_b4tp_limit_10_all
                                             , fp_b4tp_limit_0_key, fp_b4tp_limit_2_key, fp_b4tp_limit_5_key, fp_b4tp_limit_10_key
                                             )
rm(  fp_limit_0_all,      fp_limit_2_all,      fp_limit_5_all,      fp_limit_10_all
   , fp_limit_0_key,      fp_limit_2_key,      fp_limit_5_key,      fp_limit_10_key
   , fp_b4tp_limit_0_all, fp_b4tp_limit_2_all, fp_b4tp_limit_5_all, fp_b4tp_limit_10_all
   , fp_b4tp_limit_0_key, fp_b4tp_limit_2_key, fp_b4tp_limit_5_key, fp_b4tp_limit_10_key
  )

#' Rearrange columns
View(as.data.frame(names(wave_results_analysis_reshaped_fp_df)))
wave_results_analysis_reshaped_fp_df = wave_results_analysis_reshaped_fp_df[,c(  "leading_indicator_type","TFPS_cluster_min_age","TFPS_cluster_max_age","TFPS_cluster_min_descendants","LGR_p_value_threshold","parent_sub_lgr_threshold","EWS_threshold"
                                                                                ,"w2_earliest_tp_EWS_date","w2_lead_lag_days","w2_n_tp","w2_EWS_false_positive_n","w2_n_fp_before_tp","w2_n_fp_after_tp"
                                                                                ,"w3_earliest_tp_EWS_date","w3_lead_lag_days","w3_n_tp","w3_EWS_false_positive_n","w3_n_fp_before_tp","w3_n_fp_after_tp"
                                                                                ,"w4_earliest_tp_EWS_date","w4_lead_lag_days","w4_n_tp","w4_EWS_false_positive_n","w4_n_fp_before_tp","w4_n_fp_after_tp"
                                                                                ,"w5_earliest_tp_EWS_date","w5_lead_lag_days","w5_n_tp","w5_EWS_false_positive_n","w5_n_fp_before_tp","w5_n_fp_after_tp"
                                                                                ,"w6_earliest_tp_EWS_date","w6_lead_lag_days","w6_n_tp","w6_EWS_false_positive_n","w6_n_fp_before_tp","w6_n_fp_after_tp"
                                                                                ,"w7_earliest_tp_EWS_date","w7_lead_lag_days","w7_n_tp","w7_EWS_false_positive_n","w7_n_fp_before_tp","w7_n_fp_after_tp"
                                                                                ,"w8_earliest_tp_EWS_date","w8_lead_lag_days","w8_n_tp","w8_EWS_false_positive_n","w8_n_fp_before_tp","w8_n_fp_after_tp"
                                                                                ,"lead_time_total_all" , "lead_time_total_key"
                                                                                ,"false_positive_total_all" , "false_positive_total_key"
                                                                                ,"false_positive_b4tp_total_all", "false_positive_b4tp_total_key"
                                                                                ,"total_score_all", "total_score_key", "total_score_fpb4tp_all", "total_score_fpb4tp_key"
                                                                                ,"fp_limit_0_all", "fp_limit_2_all", "fp_limit_5_all", "fp_limit_10_all"
                                                                                ,"fp_limit_0_key", "fp_limit_2_key", "fp_limit_5_key", "fp_limit_10_key"
                                                                                ,"fp_b4tp_limit_0_all", "fp_b4tp_limit_2_all", "fp_b4tp_limit_5_all", "fp_b4tp_limit_10_all"
                                                                                ,"fp_b4tp_limit_0_key","fp_b4tp_limit_2_key", "fp_b4tp_limit_5_key", "fp_b4tp_limit_10_key"
                                                                              )
                                                                            ]
#' Calculate PPV and insert columns
col_spot <- which( names( wave_results_analysis_reshaped_fp_df ) == "w2_n_fp_after_tp" )
wave_results_analysis_reshaped_fp_df = data.frame(   wave_results_analysis_reshaped_fp_df[ 1 : col_spot ]
                                                  , "w2_PPV" = wave_results_analysis_reshaped_fp_df$w2_n_tp / ( wave_results_analysis_reshaped_fp_df$w2_n_tp + wave_results_analysis_reshaped_fp_df$w2_EWS_false_positive_n ) 
                                                  , "w2_PPV_fp_before_tp" = wave_results_analysis_reshaped_fp_df$w2_n_tp / ( wave_results_analysis_reshaped_fp_df$w2_n_tp + wave_results_analysis_reshaped_fp_df$w2_n_fp_before_tp ) 
                                                  , wave_results_analysis_reshaped_fp_df[ ( col_spot + 1 ) : ncol( wave_results_analysis_reshaped_fp_df ) ] 
                                                  )
col_spot <- which( names( wave_results_analysis_reshaped_fp_df ) == "w3_n_fp_after_tp" )
wave_results_analysis_reshaped_fp_df = data.frame(   wave_results_analysis_reshaped_fp_df[ 1 : col_spot ]
                                                     , "w3_PPV" = wave_results_analysis_reshaped_fp_df$w3_n_tp / ( wave_results_analysis_reshaped_fp_df$w3_n_tp + wave_results_analysis_reshaped_fp_df$w3_EWS_false_positive_n ) 
                                                     , "w3_PPV_fp_before_tp" = wave_results_analysis_reshaped_fp_df$w3_n_tp / ( wave_results_analysis_reshaped_fp_df$w3_n_tp + wave_results_analysis_reshaped_fp_df$w3_n_fp_before_tp ) 
                                                     , wave_results_analysis_reshaped_fp_df[ ( col_spot + 1 ) : ncol( wave_results_analysis_reshaped_fp_df ) ] 
                                                  )
col_spot <- which( names( wave_results_analysis_reshaped_fp_df ) == "w4_n_fp_after_tp" )
wave_results_analysis_reshaped_fp_df = data.frame(   wave_results_analysis_reshaped_fp_df[ 1 : col_spot ]
                                                     , "w4_PPV" = wave_results_analysis_reshaped_fp_df$w4_n_tp / ( wave_results_analysis_reshaped_fp_df$w4_n_tp + wave_results_analysis_reshaped_fp_df$w4_EWS_false_positive_n ) 
                                                     , "w4_PPV_fp_before_tp" = wave_results_analysis_reshaped_fp_df$w4_n_tp / ( wave_results_analysis_reshaped_fp_df$w4_n_tp + wave_results_analysis_reshaped_fp_df$w4_n_fp_before_tp ) 
                                                     , wave_results_analysis_reshaped_fp_df[ ( col_spot + 1 ) : ncol( wave_results_analysis_reshaped_fp_df ) ] 
                                                  )
col_spot <- which( names( wave_results_analysis_reshaped_fp_df ) == "w5_n_fp_after_tp" )
wave_results_analysis_reshaped_fp_df = data.frame(   wave_results_analysis_reshaped_fp_df[ 1 : col_spot ]
                                                     , "w5_PPV" = wave_results_analysis_reshaped_fp_df$w5_n_tp / ( wave_results_analysis_reshaped_fp_df$w5_n_tp + wave_results_analysis_reshaped_fp_df$w5_EWS_false_positive_n ) 
                                                     , "w5_PPV_fp_before_tp" = wave_results_analysis_reshaped_fp_df$w5_n_tp / ( wave_results_analysis_reshaped_fp_df$w5_n_tp + wave_results_analysis_reshaped_fp_df$w5_n_fp_before_tp ) 
                                                     , wave_results_analysis_reshaped_fp_df[ ( col_spot + 1 ) : ncol( wave_results_analysis_reshaped_fp_df ) ] 
                                                  )
col_spot <- which( names( wave_results_analysis_reshaped_fp_df ) == "w6_n_fp_after_tp" )
wave_results_analysis_reshaped_fp_df = data.frame(   wave_results_analysis_reshaped_fp_df[ 1 : col_spot ]
                                                     , "w6_PPV" = wave_results_analysis_reshaped_fp_df$w6_n_tp / ( wave_results_analysis_reshaped_fp_df$w6_n_tp + wave_results_analysis_reshaped_fp_df$w6_EWS_false_positive_n ) 
                                                     , "w6_PPV_fp_before_tp" = wave_results_analysis_reshaped_fp_df$w6_n_tp / ( wave_results_analysis_reshaped_fp_df$w6_n_tp + wave_results_analysis_reshaped_fp_df$w6_n_fp_before_tp ) 
                                                     , wave_results_analysis_reshaped_fp_df[ ( col_spot + 1 ) : ncol( wave_results_analysis_reshaped_fp_df ) ] 
                                                  )
col_spot <- which( names( wave_results_analysis_reshaped_fp_df ) == "w7_n_fp_after_tp" )
wave_results_analysis_reshaped_fp_df = data.frame(   wave_results_analysis_reshaped_fp_df[ 1 : col_spot ]
                                                     , "w7_PPV" = wave_results_analysis_reshaped_fp_df$w7_n_tp / ( wave_results_analysis_reshaped_fp_df$w7_n_tp + wave_results_analysis_reshaped_fp_df$w7_EWS_false_positive_n ) 
                                                     , "w7_PPV_fp_before_tp" = wave_results_analysis_reshaped_fp_df$w7_n_tp / ( wave_results_analysis_reshaped_fp_df$w7_n_tp + wave_results_analysis_reshaped_fp_df$w7_n_fp_before_tp ) 
                                                     , wave_results_analysis_reshaped_fp_df[ ( col_spot + 1 ) : ncol( wave_results_analysis_reshaped_fp_df ) ] 
                                                  )
col_spot <- which( names( wave_results_analysis_reshaped_fp_df ) == "w8_n_fp_after_tp" )
wave_results_analysis_reshaped_fp_df = data.frame(   wave_results_analysis_reshaped_fp_df[ 1 : col_spot ]
                                                     , "w8_PPV" = wave_results_analysis_reshaped_fp_df$w8_n_tp / ( wave_results_analysis_reshaped_fp_df$w8_n_tp + wave_results_analysis_reshaped_fp_df$w8_EWS_false_positive_n ) 
                                                     , "w8_PPV_fp_before_tp" = wave_results_analysis_reshaped_fp_df$w8_n_tp / ( wave_results_analysis_reshaped_fp_df$w8_n_tp + wave_results_analysis_reshaped_fp_df$w8_n_fp_before_tp ) 
                                                     , wave_results_analysis_reshaped_fp_df[ ( col_spot + 1 ) : ncol( wave_results_analysis_reshaped_fp_df ) ] 
                                                  )
#' Add mean PPVs
mean_PPV_all =  ( wave_results_analysis_reshaped_fp_df$w2_PPV + wave_results_analysis_reshaped_fp_df$w3_PPV +
                  wave_results_analysis_reshaped_fp_df$w4_PPV + wave_results_analysis_reshaped_fp_df$w5_PPV +
                  wave_results_analysis_reshaped_fp_df$w6_PPV + wave_results_analysis_reshaped_fp_df$w7_PPV +
                  wave_results_analysis_reshaped_fp_df$w8_PPV ) / 7
mean_PPV_key = (   wave_results_analysis_reshaped_fp_df$w3_PPV + wave_results_analysis_reshaped_fp_df$w4_PPV + 
                   wave_results_analysis_reshaped_fp_df$w7_PPV + wave_results_analysis_reshaped_fp_df$w8_PPV ) / 4
mean_PPV_all_fpb4tp = ( wave_results_analysis_reshaped_fp_df$w2_PPV_fp_before_tp + wave_results_analysis_reshaped_fp_df$w3_PPV_fp_before_tp +
                          wave_results_analysis_reshaped_fp_df$w4_PPV_fp_before_tp + wave_results_analysis_reshaped_fp_df$w5_PPV_fp_before_tp +
                          wave_results_analysis_reshaped_fp_df$w6_PPV_fp_before_tp + wave_results_analysis_reshaped_fp_df$w7_PPV_fp_before_tp +
                          wave_results_analysis_reshaped_fp_df$w8_PPV_fp_before_tp ) / 7
mean_PPV_key_fpb4tp = ( wave_results_analysis_reshaped_fp_df$w3_PPV_fp_before_tp + wave_results_analysis_reshaped_fp_df$w4_PPV_fp_before_tp + 
                        wave_results_analysis_reshaped_fp_df$w7_PPV_fp_before_tp + wave_results_analysis_reshaped_fp_df$w8_PPV_fp_before_tp ) / 4
wave_results_analysis_reshaped_fp_df = cbind( wave_results_analysis_reshaped_fp_df , mean_PPV_all , mean_PPV_key , mean_PPV_all_fpb4tp , mean_PPV_key_fpb4tp )

#' Save .rds and .csv
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined_2023_03_20/')
saveRDS( wave_results_analysis_reshaped_fp_df , "wave_results_analysis_reshaped_fp_df.rds" )
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/Ranking_2023_03_20/')
write.csv( wave_results_analysis_reshaped_fp_df , "wave_results_analysis_reshaped_fp_df.csv" )



#' Plot false positives against lead time for all waves and filtered by false positive limits
plot(wave_results_analysis_reshaped_fp_df$lead_time_total_all / 7 #' Mean lead lag per wave (all waves and all have true positives)
     ,wave_results_analysis_reshaped_fp_df$false_positive_total_all / 7 #' Mean false positives per wave (all waves and fp before and after 1st tp)
     ,ylim=c(0,10)
     ,xlab="Mean earliest EWS lead (-ve) or lag (+ve) days per wave"
     ,ylab="Mean EWS false positives per wave")
#' Add points for different false positive filters
sub_fp_limit_10_all = subset( wave_results_analysis_reshaped_fp_df , wave_results_analysis_reshaped_fp_df$fp_limit_10_all == "y" )
sub_fp_limit_5_all  = subset( wave_results_analysis_reshaped_fp_df , wave_results_analysis_reshaped_fp_df$fp_limit_5_all  == "y" )
sub_fp_limit_2_all  = subset( wave_results_analysis_reshaped_fp_df , wave_results_analysis_reshaped_fp_df$fp_limit_2_all  == "y" )
sub_fp_limit_0_all  = subset( wave_results_analysis_reshaped_fp_df , wave_results_analysis_reshaped_fp_df$fp_limit_0_all  == "y" )
#' false positives <= 10 
points( sub_fp_limit_10_all$lead_time_total_all / 7 , sub_fp_limit_10_all$false_positive_total_all / 7 , col="blue"   , pch=16 )
#' false positives <= 5 
points( sub_fp_limit_5_all$lead_time_total_all   / 7 , sub_fp_limit_5_all$false_positive_total_all / 7 , col="grey"   , pch=16 )
#' false positives <= 2 
points( sub_fp_limit_2_all$lead_time_total_all   / 7 , sub_fp_limit_2_all$false_positive_total_all / 7 , col="orange" , pch=16 )
#' false positives = 0
points( sub_fp_limit_0_all$lead_time_total_all   / 7 , sub_fp_limit_0_all$false_positive_total_all / 7 , col="red"    , pch=16 )
legend("topright"
       ,legend=c("False positive EWS All","False positive EWS <=10","False positive EWS <=5","False positive EWS <=2","False positive EWS =0")
       ,col=c("black","blue","grey","orange","red")
       ,pch=16)

#' Histograms to show shift in lead time with different limits on number of false positives
hist(  sub_fp_limit_10_all$lead_time_total_all / 7 #' Mean lead lag per wave (all waves and all have true positives)
       , breaks = 20
       , col=rgb(0,0,1,0.5) # red,green,blue and transparency
       , xlim = c(-10,40)
       , ylim = c(0,0.3)
       , freq = FALSE
       , xlab = "EWS lead (-ve) or lag (+ve) days"
       #, main = paste("Impact of parent/sub-cluster replacement on EWS lead/lag for Omicron BA.1 wave"))
       #, main = paste("No parent/sub-cluster replacement"))
       , main = paste("Distribution of lead times changes with number of false positives allowed in filter"))
hist(sub_fp_limit_5_all$lead_time_total_all / 7, breaks = 10, col=rgb(0,1,0,0.5), add=T, freq = FALSE)
hist(sub_fp_limit_2_all$lead_time_total_all / 7, breaks = 10, col=rgb(1,0,0,0.5), add=T, freq = FALSE)
hist(sub_fp_limit_0_all$lead_time_total_all / 7, breaks = 40, col=rgb(1,1,0,0.5), add=T, freq = FALSE)
#hist(na.omit(merged_df_not_999_100$w7_lead_lag_days), breaks = 20, col=rgb(1,1,0,0.5), add=T, freq = FALSE)
legend("topright"
       ,legend=c("False positive EWS <=10","False positive EWS <=5","False positive EWS <=2","False positive EWS =0")
       ,col=c(rgb(0,0,1,0.5),rgb(0,1,0,0.5),rgb(1,0,0,0.5),rgb(1,1,0,0.5))
       ,lty=1,lwd=5)
#' Histogram seems to show that as tighten the filter on false positives the lead time reduces. This is counter intuitive. 
#' Possibly due to the different leading indicators included in the different sets. 
#' Therefore, need to analyse based on individual leading indicators.

#' Plot means of means so get a single point instead of a cloud of all parameter sets
plot(  mean( wave_results_analysis_reshaped_fp_df$lead_time_total_all / 7 ) #' Mean lead lag per wave (all waves and all have true positives)
      ,mean( wave_results_analysis_reshaped_fp_df$false_positive_total_all / 7 ) #' Mean false positives per wave (all waves and fp before and after 1st tp)
      ,xlim=c(-10,30)
      ,ylim=c(0,10)
      ,xlab="Mean earliest EWS lead (-ve) or lag (+ve) days per wave"
      ,ylab="Mean EWS false positives per wave"
      ,pch=16, cex=2)
#' false positives <= 10 
points(   mean( sub_fp_limit_10_all$lead_time_total_all / 7 ) #' Mean lead lag per wave (all waves and all have true positives)
        , mean( sub_fp_limit_10_all$false_positive_total_all / 7 ) #' Mean false positives per wave (all waves and fp before and after 1st tp)
        , col="blue"
        , pch =16 , cex=2)
#' false positives <= 5 
points(     mean( sub_fp_limit_5_all$lead_time_total_all / 7 ) #' Mean lead lag per wave (all waves and all have true positives)
          , mean( sub_fp_limit_5_all$false_positive_total_all / 7 ) #' Mean false positives per wave (all waves and fp before and after 1st tp)
          , col="grey"
          , pch =16 , cex=2)
#' false positives <= 2 
points(     mean( sub_fp_limit_2_all$lead_time_total_all / 7 ) #' Mean lead lag per wave (all waves and all have true positives)
          , mean( sub_fp_limit_2_all$false_positive_total_all / 7 ) #' Mean false positives per wave (all waves and fp before and after 1st tp)
          , col="orange"
          , pch =16 , cex=2 )
#' false positives = 0
points(   mean( sub_fp_limit_0_all$lead_time_total_all / 7 ) #' Mean lead lag per wave (all waves and all have true positives)
        , mean( sub_fp_limit_0_all$false_positive_total_all / 7 ) #' Mean false positives per wave (all waves and fp before and after 1st tp)
        , col="red"
        , pch =16 , cex=2 )

#hist(na.omit(merged_df_not_999_100$w7_lead_lag_days), breaks = 20, col=rgb(1,1,0,0.5), add=T, freq = FALSE)
legend("topright"
       ,legend=c("Any number of false positive EWS","False positive EWS <=10","False positive EWS <=5","False positive EWS <=2","False positive EWS =0")
       ,col=c("black","blue","grey","orange","red")
       ,pch=16)

#' looks like the shift is actually in the other direction to that expected (bigger lead time with more stringent limit on false positives)
#' Perhaps it is due to the tendency of leading indicators to have better lead times or more/less false positives
#' As reduce the number limit on false positives pango_lineage leading indicator becomes a higher proportion and it has better lead times
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/leading_indicator_proportions_2023_03_20")
write.csv( round( table( sub_fp_limit_10_all$leading_indicator_type ) / sum( table( sub_fp_limit_10_all$leading_indicator_type ) ) * 100 , 1 ) , "fp_limit_10_li_proportion.csv" )
write.csv( round( table( sub_fp_limit_5_all$leading_indicator_type )  / sum( table( sub_fp_limit_5_all$leading_indicator_type  ) ) * 100 , 1 ) , "fp_limit_5_li_proportion.csv"  )
write.csv( round( table( sub_fp_limit_2_all$leading_indicator_type )  / sum( table( sub_fp_limit_2_all$leading_indicator_type  ) ) * 100 , 1 ) , "fp_limit_2_li_proportion.csv"  )
write.csv( round( table( sub_fp_limit_0_all$leading_indicator_type )  / sum( table( sub_fp_limit_0_all$leading_indicator_type  ) ) * 100 , 1 ) , "fp_limit_0_li_proportion.csv"  ) 


#' Repeat but only for 'key' waves driven by genomic variants: Alpha, Delta 1, Omicron BA.1 and Omicron BA.2 
#' (Note that B.1.177 excluded because limited data in run up to the wave)
#' Note that this still uses the requirement that there must be a True Positive EWS for all waves. 
#' It is just the ranking that uses these 'key' waves and also the calculation of the mean earliest EWS lead time and 
#' mean EWS false positive

plot( wave_results_analysis_reshaped_fp_df$lead_time_total_key / 4 #' Mean lead lag per wave (all waves and all have true positives)
     ,wave_results_analysis_reshaped_fp_df$false_positive_total_key / 4 #' Mean false positives per wave (all waves and fp before and after 1st tp)
     #,ylim=c(0,10)
     ,xlab="Mean earliest EWS lead (-ve) or lag (+ve) days per wave"
     ,ylab="Mean EWS false positives per wave")
#' Add points for different false positive filters
sub_fp_limit_10_key = subset( wave_results_analysis_reshaped_fp_df , wave_results_analysis_reshaped_fp_df$fp_limit_10_key == "y" )
sub_fp_limit_5_key  = subset( wave_results_analysis_reshaped_fp_df , wave_results_analysis_reshaped_fp_df$fp_limit_5_key == "y" )
sub_fp_limit_2_key  = subset( wave_results_analysis_reshaped_fp_df , wave_results_analysis_reshaped_fp_df$fp_limit_2_key == "y" )
sub_fp_limit_0_key  = subset( wave_results_analysis_reshaped_fp_df , wave_results_analysis_reshaped_fp_df$fp_limit_0_key == "y" )
#' false positives <= 10 
points(   sub_fp_limit_10_key$lead_time_total_key / 4 #' Mean lead lag per wave (all waves and all have true positives)
        , sub_fp_limit_10_key$false_positive_total_key / 4 #' Mean false positives per wave (all waves and fp before and after 1st tp)
        , col="blue"
        , pch =16)
#' false positives <= 5 
points(   sub_fp_limit_5_key$lead_time_total_key / 4 #' Mean lead lag per wave (all waves and all have true positives)
        , sub_fp_limit_5_key$false_positive_total_key / 4 #' Mean false positives per wave (all waves and fp before and after 1st tp)
        , col="grey"
        , pch =16)
#' false positives <= 2 
points(   sub_fp_limit_2_key$lead_time_total_key / 4 #' Mean lead lag per wave (all waves and all have true positives)
        , sub_fp_limit_2_key$false_positive_total_key / 4 #' Mean false positives per wave (all waves and fp before and after 1st tp)
        , col="orange"
        , pch =16)
#' false positives = 0
points(   sub_fp_limit_0_key$lead_time_total_key / 4 #' Mean lead lag per wave (all waves and all have true positives)
        , sub_fp_limit_0_key$false_positive_total_key / 4 #' Mean false positives per wave (all waves and fp before and after 1st tp)
        , col="red"
        , pch =16)
legend("topright"
       ,legend=c("Any number of false positive EWS","False positive EWS <=10","False positive EWS <=5","False positive EWS <=2","False positive EWS =0")
       ,col=c("black","blue","grey","orange","red")
       ,pch=16)

#' Histograms to show shift in lead time with different limits on number of false positives
hist(  sub_fp_limit_10_key$lead_time_total_key / 4 #' Mean lead lag per wave (all waves and all have true positives)
       , breaks = 20
       , col=rgb(0,0,0,0.5) # red,green,blue and transparency
       #, xlim = c(-10,40)
       , ylim = c(0,0.3)
       , freq = FALSE
       , xlab = "EWS lead (-ve) or lag (+ve) days"
       #, main = paste("Impact of parent/sub-cluster replacement on EWS lead/lag for Omicron BA.1 wave"))
       #, main = paste("No parent/sub-cluster replacement"))
       , main = paste("Distribution of lead times changes with number of false positives allowed in filter"))
hist(sub_fp_limit_5_key$lead_time_total_key / 4, breaks = 20, add=T, freq = FALSE, col=rgb(0,1,0,0.5) )
hist(sub_fp_limit_2_key$lead_time_total_key / 4, breaks = 20, add=T, freq = FALSE, col=rgb(1,0.5,0,0.5))
hist(sub_fp_limit_0_key$lead_time_total_key / 4, breaks = 20, add=T, freq = FALSE, col=rgb(1,0,0,0.5))
legend("topright"
       ,legend=c("False positive EWS <=10","False positive EWS <=5","False positive EWS <=2","False positive EWS =0")
       ,col=c(rgb(0,0,0,0.5),rgb(0,1,0,0.5),rgb(1,0.5,0,0.5),rgb(1,0,0,0.5))
       ,lty=1,lwd=5)

#' Difficult to determine shift. Distribution looks bi-modal. 
#' But definite shift to bigger lead time (smaller lag) when false positives limited to zero.
#' Again, this is counter intuitive.
#' Looks like it is due to pango leading indicator having less false positives (which makes sense) and having a good lead time
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/leading_indicator_proportions_2023_03_20")
write.csv( round( table( wave_results_analysis_reshaped_fp_df$leading_indicator_type ) / sum( table( wave_results_analysis_reshaped_fp_df$leading_indicator_type ) ) * 100 , 1 ) , "fp_limit_all_li_proportion.csv" )
write.csv( round( table( sub_fp_limit_10_key$leading_indicator_type ) / sum( table( sub_fp_limit_10_key$leading_indicator_type ) ) * 100 , 1 ) , "fp_limit_10_key_li_proportion.csv" )
write.csv( round( table( sub_fp_limit_5_key$leading_indicator_type )  / sum( table( sub_fp_limit_5_key$leading_indicator_type ) )  * 100 , 1 ) , "fp_limit_5_key_li_proportion.csv" )
write.csv( round( table( sub_fp_limit_2_key$leading_indicator_type )  / sum( table( sub_fp_limit_2_key$leading_indicator_type ) )  * 100 , 1 ) , "fp_limit_2_key_li_proportion.csv" )
write.csv( round( table( sub_fp_limit_0_key$leading_indicator_type )  / sum( table( sub_fp_limit_0_key$leading_indicator_type ) )  * 100 , 1 ) , "fp_limit_0_key_li_proportion.csv" )


#' plot scores: All waves vs Key waves
par( mfrow = c( 1 , 1 ) )
plot( wave_results_analysis_reshaped_fp_df$lead_time_total_all / 7
      ,wave_results_analysis_reshaped_fp_df$false_positive_total_all / 7
      ,ylim=c(0,12)
      ,xlim=c(-20,40)
      #,col="red"
      , pch=16
      ,xlab="Mean earliest EWS lead/lag days per wave"
      ,ylab="Mean EWS false positives per wave"
      ,cex=1.2)
points( wave_results_analysis_reshaped_fp_df$lead_time_total_key / 4
       ,wave_results_analysis_reshaped_fp_df$false_positive_total_key / 4 
       ,col="red",pch=16)
#,ylim=c(0,10)
#,xlab="average earliest EWS lead/lag days per wave"
#,ylab="average EWS false positives per wave")
abline(h=0)
abline(v=0)
legend( "topright", legend = c("All UK infection waves", "Alpha, Delta 1, Omicron BA.1, Omicron BA.2"), col=c("black","red") , pch=16)
legend( "topright", legend = c("All leading indicators - All UK infection waves", "All leading indicators - Alpha, Delta 1, Omicron BA.1, Omicron BA.2", "Pango lineage LGR max - All UK infection waves","Pango lineage LGR max - Alpha, Delta 1, Omicron BA.1, Omicron BA.2"), col=c("black","red","blue","lightblue") , pch=16)
subset_pango_max = subset( wave_results_analysis_reshaped_fp_df, wave_results_analysis_reshaped_fp_df$leading_indicator_type == "tfps_pango_lineage_max_lgr")
points(    subset_pango_max$lead_time_total_all / 7
         , subset_pango_max$false_positive_total_all / 7
         , col="blue")
points(    subset_pango_max$lead_time_total_key / 4
         , subset_pango_max$false_positive_total_key / 4
         , col="lightblue"
         , pch=16)

#' Plot means of means
par( mfrow = c( 1 , 1 ) )
plot(    mean( wave_results_analysis_reshaped_fp_df$lead_time_total_all / 7 )
       , mean( wave_results_analysis_reshaped_fp_df$false_positive_total_all / 7 )
       , ylim=c(0,10) , xlim=c(-10,50)
       #,col="red"
       , xlab="Mean earliest EWS lead/lag days per wave", ylab="Mean EWS false positives per wave"
       , pch=16,cex=2)
points(    mean( wave_results_analysis_reshaped_fp_df$lead_time_total_key / 4 )
         , mean( wave_results_analysis_reshaped_fp_df$false_positive_total_key / 4 )
         , col="red", pch=16, cex=2)
points(    median( wave_results_analysis_reshaped_fp_df$lead_time_total_key / 4 )
         , median( wave_results_analysis_reshaped_fp_df$false_positive_total_key / 4 )
         , col="pink", pch=16, cex=2)

abline(h=0); abline(v=0)
legend( "topright", legend = c("All UK infection waves", "Alpha, Delta 1, Omicron BA.1, Omicron BA.2"), col=c("black","red") , pch=1)
legend( "topright", legend = c("All leading indicators - All UK infection waves", "All leading indicators - Alpha, Delta 1, Omicron BA.1, Omicron BA.2", "Pango lineage LGR max - All UK infection waves","Pango lineage LGR max - Alpha, Delta 1, Omicron BA.1, Omicron BA.2"), col=c("black","red","blue","lightblue") , pch=16)
subset_pango_max = subset( wave_results_analysis_reshaped_fp_df, wave_results_analysis_reshaped_fp_df$leading_indicator_type == "tfps_pango_lineage_max_lgr")
points(    mean( subset_pango_max$lead_time_total_all / 7 )
         , mean( subset_pango_max$false_positive_total_all / 7 )
         , col="blue" , pch=16 , cex=2)
points(  median( subset_pango_max$lead_time_total_all / 7 )
         , median( subset_pango_max$false_positive_total_all / 7 )
         , col="purple" , pch=16 , cex=2)
points(  mean( subset_pango_max$lead_time_total_key / 4 )
         , mean( subset_pango_max$false_positive_total_key / 4 )
         , col="lightblue" , pch=16 , cex =2)
points(  median( subset_pango_max$lead_time_total_key / 4 )
         , median( subset_pango_max$false_positive_total_key / 4 )
         , col="lightgreen" , pch=16 , cex =2)


#' Plot mean EWS false positives per wave against mean earliest EWS lead (-ve) or lag (+ve) per wave by leading indicator
#' loop through all leading indicators in dataframe and create a plot for each leading indicator
for (li in unique( wave_results_analysis_reshaped_fp_df$leading_indicator_type )){ #[-1]){
  #' Reset plot values
  x_all = NULL ; x_key = NULL ; m_all = NULL ; m_key = NULL
  #' Subset by required leading indicator (by changing the index value)
  li_filtered = subset( wave_results_analysis_reshaped_fp_df 
                        , wave_results_analysis_reshaped_fp_df$leading_indicator_type == li )
  x_all = li_filtered$lead_time_total_all/7 ; y_all = li_filtered$false_positive_total_all/7
  x_key = li_filtered$lead_time_total_key/4 ; y_key = li_filtered$false_positive_total_key/4
  #' Format data for fitting
  df_all = data.frame( x_all , y_all )
  df_all = subset( df_all , !is.na(df_all$x_all) & !is.na(df_all$y_all) )
  dt_all = data.table::setorder( df_all , "x_all" )
  x_all = dt_all$x_all
  y_all = dt_all$y_all
  df_key = data.frame( x_key , y_key )
  df_key = subset( df_key , !is.na( df_key$x_key ) & !is.na( df_key$y_key ) )
  dt_key = data.table::setorder( df_key , "x_key" )
  x_key = dt_key$x_key
  y_key = dt_key$y_key
  #' Linear fit to scatter plots
  fit_all <- lm( y_all ~ x_all ) 
  fit_key <- lm( y_key ~ x_key )
  #' GAM fit
  library(mgcv)
  GAM_smooth_function ="ps" ;  deg_free_k_all = min( 11 , nrow(df_all) )-1 ; deg_free_k_key = min( 11 , nrow(df_key) )-1
  m_all <- try( mgcv::gam( y_all ~ s( x_all, bs=GAM_smooth_function, k = deg_free_k_all), data = dt_all ) , silent = TRUE )
  #if ('try-error' %in% class(m_all)) next
  m_key <- try( mgcv::gam( y_key ~ s( x_key, bs=GAM_smooth_function, k = deg_free_k_key), data = dt_key ) , silent = TRUE )
  #if ('try-error' %in% class(m_key)) next
  #' Plot data
  setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/plots_2023_03_21")
  png(paste0(li,"_GAM_k_",deg_free_k_all,".png"),width=9.64,height=6.35,units="in",res=400) #' save to file
  plot( x_all , y_all , ylim = c( 0 , 20 ) , xlim = c( -20 , 50 ), main = li
        ,xlab="Mean earliest EWS lead (-ve) or lag (+ve) days per wave" , ylab="Mean EWS false positives per wave"
        , col = "grey"
        , xaxs = "i" , yaxs = "i")
  points( x_key , y_key , col = "green" )
  abline(h=0)
  abline(v=0)
  legend( "topright"
          , legend = c("Data: All UK infection waves"
                       , paste0("GAM(ps,k=",deg_free_k_all,"): All UK infection waves")
                       , "Linear: All UK infection waves"
                       , "Data: Alpha, Delta 1, Omicron BA.1, Omicron BA.2"
                       , paste0("GAM(ps,k=",deg_free_k_key,"): Alpha, Delta 1, Omicron BA.1, Omicron BA.2") 
                       , "Linear: Alpha, Delta 1, Omicron BA.1, Omicron BA.2" )
          , col = c("grey","darkgrey","darkgrey","lightgreen","darkgreen","darkgreen") 
          , pch = c( 1,NA,NA, 1,NA,NA)
          , lty = c(NA, 1, 2, NA, 1, 2)
          , lwd = c( 1, 2, 2, 1, 2, 2))
  #' Plot linear models
  lines( x_all , predict( fit_all ) , col = "darkgrey"  , lwd= 2 , lty =2 )  #lines( x_axis_all , predict( fit_all , data.frame( x = x_axis_all ) ) , col = "black" )
  lines( x_key , predict( fit_key ) , col = "darkgreen" , lwd= 2 , lty =2 )
  #' Plot Generalised Additive models
  try( lines( m_all$model[[2]] , m_all$fitted.values , typ="l" , col = "darkgrey" , lwd=2 ), silent = TRUE )
  try( lines( m_key$model[[2]] , m_key$fitted.values , typ="l" , col = "darkgreen" , lwd=2 ), silent = TRUE )
  #' Turn off saving plot to file
  dev.off()
}

#' Plot mean EWS false positives per wave against mean earliest EWS lead (-ve) or lag (+ve) per wave by leading indicator
#' loop through all leading indicators in dataframe and add the linear fit for each leading indicator to a single plot
plot( 0 
      , 0 
      , xlim=c(-20,50), ylim=c(0,10)
      , xaxs = "i" , yaxs = "i"
      , xlab="mean earliest EWS lead (-ve) or lag (+ve) days per wave" , ylab="mean EWS false positives per wave"
      , col = "darkgrey"  , lwd= 2 , lty =1 )  #lines( x_axis_all , predict( fit_all , data.frame( x = x_axis_all ) ) , col = "black" )
abline(v=0)
for (li in unique( wave_results_analysis_reshaped_fp_df$leading_indicator_type )){ #[-1]){
  #' Reset plot values
  x_all = NULL ; x_key = NULL ; m_all = NULL ; m_key = NULL
  #' Subset by required leading indicator (by changing the index value)
  li_filtered = subset( wave_results_analysis_reshaped_fp_df 
                        , wave_results_analysis_reshaped_fp_df$leading_indicator_type == li )
  x_all = li_filtered$lead_time_total_all/7 ; y_all = li_filtered$false_positive_total_all/7
  x_key = li_filtered$lead_time_total_key/4 ; y_key = li_filtered$false_positive_total_key/4
  #' Format data for fitting
  df_all = data.frame( x_all , y_all )
  df_all = subset( df_all , !is.na(df_all$x_all) & !is.na(df_all$y_all) )
  dt_all = data.table::setorder( df_all , "x_all" )
  x_all = dt_all$x_all
  y_all = dt_all$y_all
  df_key = data.frame( x_key , y_key )
  df_key = subset( df_key , !is.na( df_key$x_key ) & !is.na( df_key$y_key ) )
  dt_key = data.table::setorder( df_key , "x_key" )
  x_key = dt_key$x_key
  y_key = dt_key$y_key
  #' Linear fit to scatter plots
  fit_all <- lm( y_all ~ x_all ) 
  fit_key <- lm( y_key ~ x_key )
  #' Plot data
  #' Plot linear models
  lines( x_all 
        , predict( fit_all ) 
        , xlim=c(-20,50), ylim=c(0,10)
        , xaxs = "i" , yaxs = "i"
        , xlab="mean earliest EWS lead (-ve) or lag (+ve) days per wave" , ylab="mean EWS false positives per wave"
        , col = "blue"  , lwd= 2 , lty =1 )  #lines( x_axis_all , predict( fit_all , data.frame( x = x_axis_all ) ) , col = "black" )
  lines( x_key 
          , predict( fit_key ) 
          , xaxs = "i" , yaxs = "i"
          , col = "red"  , lwd= 2 , lty =1 )  #lines( x_axis_all , predict( fit_all , data.frame( x = x_axis_all ) ) , col = "black" )
  legend( "topright"
        , legend = c(  "Linear: All UK infection waves"
                     , "Linear: Alpha, Delta 1, Omicron BA.1, Omicron BA.2" )
        , col = c("blue","red") 
        #, lty = c(2, 2)
        , lwd = c( 2, 2))
}

################################################
#' Repeat analysis above but only using number of false positive EWS BEFORE the earliest true positive instead of all
################################################
#' Rating the ability of different parameter sets in generating EWS
#' using total number of false positives EWS BEFORE the earliest true positive instead of all

plot(    wave_results_analysis_reshaped_fp_df$lead_time_total_all / 7 #' Only divide by 7 if filtered for parameter sets with a true positive for all 7 waves
       , wave_results_analysis_reshaped_fp_df$false_positive_b4tp_total_all / 7 
       , xlim=c(-20,40 ), ylim=c(0,12)
       , xlab="Mean earliest EWS lead (-ve) or lag (+ve) days per wave"
       , ylab="Mean EWS false positives (before earliest true positive) per wave"
       ,pch=16)
points( wave_results_analysis_reshaped_fp_df$lead_time_total_key / 4 #' Mean lead or lag for key waves
       ,wave_results_analysis_reshaped_fp_df$false_positive_b4tp_total_key / 4 #' Mean false positives before true positives for key waves
       ,col="red",pch=16)
abline(h=0)
abline(v=0)
legend( "topright", legend = c("All leading indicators - All UK infection waves", "All leading indicators - Alpha, Delta 1, Omicron BA.1, Omicron BA.2"), col=c("black","red") , pch=16)
legend( "topright", legend = c("All leading indicators - All UK infection waves", "All leading indicators - Alpha, Delta 1, Omicron BA.1, Omicron BA.2", "Pango lineage LGR max - All UK infection waves","Pango lineage LGR max - Alpha, Delta 1, Omicron BA.1, Omicron BA.2"), col=c("black","red","blue","lightblue") , pch=16)
subset_pango_max = subset( wave_results_analysis_reshaped_fp_df , wave_results_analysis_reshaped_fp_df$leading_indicator_type == "tfps_pango_lineage_max_lgr")
points(    subset_pango_max$lead_time_total_all / 7
         , subset_pango_max$false_positive_b4tp_total_all / 7
         , col="blue")
points(  subset_pango_max$lead_time_total_key / 4
         , subset_pango_max$false_positive_b4tp_total_key / 4
         , col="lightblue"
         , pch=16)

#' Histogram of false positive EWS per wave for parameter sets producing at least one true positive for each of 7 waves
hist( as.numeric( as.vector( wave_results_analysis_reshaped_fp_df$false_positive_b4tp_total_all / 7 ) ) , breaks = 20)
hist( as.numeric( as.vector( wave_results_analysis_reshaped_fp_df$lead_time_total_all / 7 ) ) , breaks = 20)

#' Plot mean EWS false positives (BEFORE earliest tp) per wave against mean earliest EWS lead (-ve) or lag (+ve) per wave by leading indicator
#' loop through all leading indicators in datafarame and create a plot for each leading indicator
for (li in unique( wave_results_analysis_reshaped_fp_df$leading_indicator_type )){
  #' Reset plot values
  x_all = NULL ; x_key = NULL ; m_all = NULL ; m_key = NULL
  #' Subset by required leading indicator (by changing the index value)
  li_filtered = subset( wave_results_analysis_reshaped_fp_df 
                        , wave_results_analysis_reshaped_fp_df$leading_indicator_type == li )
  x_all = li_filtered$lead_time_total_all / 7 ; y_all = li_filtered$false_positive_b4tp_total_all / 7
  x_key = li_filtered$lead_time_total_key / 4 ; y_key = li_filtered$false_positive_b4tp_total_key / 4
  #' Format data for fitting
  df_all = data.frame( x_all , y_all )
  df_all = subset( df_all , !is.na(df_all$x_all) & !is.na(df_all$y_all) )
  dt_all = data.table::setorder( df_all , "x_all" )
  x_all = dt_all$x_all
  y_all = dt_all$y_all
  df_key = data.frame( x_key , y_key )
  df_key = subset( df_key , !is.na( df_key$x_key ) & !is.na( df_key$y_key ) )
  dt_key = data.table::setorder( df_key , "x_key" )
  x_key = dt_key$x_key
  y_key = dt_key$y_key
  #' Linear fit to scatter plots
  fit_all <- lm( y_all ~ x_all ) 
  fit_key <- lm( y_key ~ x_key )
  #' GAM fit
  library(mgcv)
  GAM_smooth_function ="ps" ;  deg_free_k_all = min( 11 , nrow(df_all) )-1 ; deg_free_k_key = min( 11 , nrow(df_key) )-1
  m_all <- try( mgcv::gam( y_all ~ s( x_all, bs=GAM_smooth_function, k = deg_free_k_all), data = dt_all ) , silent = TRUE )
  #if ('try-error' %in% class(m_all)) next
  m_key <- try( mgcv::gam( y_key ~ s( x_key, bs=GAM_smooth_function, k = deg_free_k_key), data = dt_key ) , silent = TRUE )
  #if ('try-error' %in% class(m_key)) next
  #' Plot data
  setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/plots_2023_03_21/")
  png(paste0(li,"_GAM_k_",deg_free_k_all,"fp_before_tp.png"),width=9.64,height=6.35,units="in",res=400) #' save to file
  plot(   x_all , y_all , ylim = c( 0 , 20 ) , xlim = c( -20 , 50 ), main = li
        , xlab="Mean earliest EWS lead (-ve) or lag (+ve) days per wave" , ylab="Mean EWS false positives (before earliest true positive) per wave"
        , col = "grey"
        , xaxs = "i" , yaxs = "i")
  points( x_key , y_key , col = "green" )
  abline(h=0)
  abline(v=0)
  legend( "topright"
          , legend = c("Data: All UK infection waves"
                       , paste0("GAM(ps,k=",deg_free_k_all,"): All UK infection waves")
                       , "Linear: All UK infection waves"
                       , "Data: Alpha, Delta 1, Omicron BA.1, Omicron BA.2"
                       , paste0("GAM(ps,k=",deg_free_k_key,"): Alpha, Delta 1, Omicron BA.1, Omicron BA.2") 
                       , "Linear: Alpha, Delta 1, Omicron BA.1, Omicron BA.2" )
          , col = c("grey","darkgrey","darkgrey","lightgreen","darkgreen","darkgreen") 
          , pch = c( 1,NA,NA, 1,NA,NA)
          , lty = c(NA, 1, 2, NA, 1, 2)
          , lwd = c( 1, 2, 2, 1, 2, 2))
  #' Plot linear models
  lines( x_all , predict( fit_all ) , col = "darkgrey" , lwd= 2 , lty =2 )  #lines( x_axis_all , predict( fit_all , data.frame( x = x_axis_all ) ) , col = "black" )
  lines( x_key , predict( fit_key ) , col = "darkgreen" , lwd= 2 , lty =2 )
  #' Plot Generalised Additive models
  try( lines( m_all$model[[2]] , m_all$fitted.values , typ="l" , col = "darkgrey" , lwd=2 ), silent = TRUE )
  try( lines( m_key$model[[2]] , m_key$fitted.values , typ="l" , col = "darkgreen" , lwd=2 ), silent = TRUE )
  #' Turn off saving plot to file
  dev.off()
  #invisible(readline(prompt="Press [enter] to continue"))
}



#############
#' Plot top ranked leading indicator, hospitalisations, cases, EWS threshold, true positives and false positives 
#############
#' Load case/hospitalisation data
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs'
#filename <- "data_2022-May-09 - UK Covid-19 cases by sample date.csv"
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
#dat_type <- "cases"
dat_type <- "hospitalisations"
#' data_load() @ C:\Users\kdrake\GitHub\Early_Warning_Signal
library(data.table)
dat_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)
par( mfrow = c( 3 , 1 ) )
plot( dat_df$date
      , dat_df$cases
      , xlab = "Date"
      , ylab = "UK Covid-19 hospitalisations"
      , cex.lab = 1.25
      , cex.axis = 1.25
      , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )
      , ylim = c(0, 5000)
      #, xaxt = "n"
      , xaxs="i"
      , yaxs="i" #' Makes plot box tight to values
      , typ="l"
      , lwd = 1.5
      , col = "black" #data_colour[2]
)
#' Wave start dates (average where a wave has more than one start date as 
#' produced from the 12 models incorporating 'low resolution wave filter' - see
#' 21 June 2022 update slides)
wave_start_dates <- c(            #' 1 Wuhan
  as.Date("2020-08-19") #' 2 B.1.177
  ,as.Date("2020-11-29") #' 3 Alpha
  ,as.Date("2021-05-11") #' 4 Delta
  ,as.Date("2021-08-03") #' 5 Delta
  ,as.Date("2021-09-27") #' 6 Delta
  ,as.Date("2021-11-26") #' 7 Omicron
  ,as.Date("2022-02-21") #' 8 Omicron
  #' 9 Omicron
)

#' Read in case data?
#lines()
#' Read in dataframe containing information on early warning signals (as calculated above) 
#' for each set of scan variables (min age, max age, min descendants) and filter variable (LGR p-value) 
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined_2023_03_20/")
ews_results_df = readRDS( "ews_results_df.rds" )
parameter_set_ews_results_df = subset( ews_results_df , ews_results_df$TFPS_cluster_min_age == "07" & #"07" 14" "28"
                                                    ews_results_df$TFPS_cluster_max_age == "56" & #"56" 84"
                                                    ews_results_df$TFPS_cluster_min_descendants == "20" & #"050" "100" "proportional"
                                                    ews_results_df$LGR_p_value_threshold == 0.01 & #0.05 #0.01 #1000
                                                    ews_results_df$parent_sub_lgr_threshold == 0.85 & #seq(0.60,1.00,0.05) and 999
                                                    ews_results_df$leading_indicator_type == "tfps_pango_lineage_max_lgr" &
                                                    ews_results_df$EWS_threshold == 0.00 #' seq(0.00,5.00,0.05)
                                  )
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined_2023_03_20/")
wave_results_df = readRDS( "wave_results_df.rds" )
View(wave_results_df[1:10,])
parameter_set_wave_results_df = subset( wave_results_df , wave_results_df$TFPS_cluster_min_age == "07" & #"07" 14" "28"
                                     wave_results_df$TFPS_cluster_max_age == "56" & #"56" 84"
                                     wave_results_df$TFPS_cluster_min_descendants == "20" & #"050" "100" "proportional"
                                     wave_results_df$LGR_p_value_threshold == 10000 & #0.05 #0.01
                                     wave_results_df$parent_sub_LGR_threshold == 0.90 & #seq(0.60,1.00,0.05) and 999
                                     wave_results_df$leading_indicator_type == "tfps_pango_lineage_max_lgr" &
                                     wave_results_df$EWS_threshold == 0.45 #' seq(0.00,5.00,0.05)
                                      )
true_positive = subset( parameter_set_wave_results_df , parameter_set_wave_results_df$positive_EWS == TRUE )
false_positive = subset( parameter_set_wave_results_df , parameter_set_wave_results_df$positive_EWS == FALSE ) 
parameter_set_wave_results_analysis_reshaped_df = subset( wave_results_analysis_reshaped_df , wave_results_analysis_reshaped_df$TFPS_cluster_min_age == "07" & #"07" 14" "28"
                                                                                              wave_results_analysis_reshaped_df$TFPS_cluster_max_age == "56" & #"56" 84"
                                                                                              wave_results_analysis_reshaped_df$TFPS_cluster_min_descendants == "20" & #"050" "100" "proportional"
                                                                                              wave_results_analysis_reshaped_df$LGR_p_value_threshold == 10000 & #0.05 #0.01
                                                                                              wave_results_analysis_reshaped_df$parent_sub_lgr_threshold == 0.90 & #seq(0.60,1.00,0.05) and 999
                                                                                              wave_results_analysis_reshaped_df$leading_indicator_type == "tfps_pango_lineage_max_lgr" &
                                                                                              wave_results_analysis_reshaped_df$EWS_threshold == 0.45 #' seq(0.00,5.00,0.05)
                                                          )
true_positive_earliest = c(   as.Date( parameter_set_wave_results_analysis_reshaped_df$w2_earliest_tp_EWS_date , origin = "1970-01-01" )
                            , as.Date( parameter_set_wave_results_analysis_reshaped_df$w3_earliest_tp_EWS_date , origin = "1970-01-01" )
                            , as.Date( parameter_set_wave_results_analysis_reshaped_df$w4_earliest_tp_EWS_date , origin = "1970-01-01" )
                            , as.Date( parameter_set_wave_results_analysis_reshaped_df$w5_earliest_tp_EWS_date , origin = "1970-01-01" )
                            , as.Date( parameter_set_wave_results_analysis_reshaped_df$w6_earliest_tp_EWS_date , origin = "1970-01-01" )
                            , as.Date( parameter_set_wave_results_analysis_reshaped_df$w7_earliest_tp_EWS_date , origin = "1970-01-01" )
                            , as.Date( parameter_set_wave_results_analysis_reshaped_df$w8_earliest_tp_EWS_date , origin = "1970-01-01" )
                          )
#' Plot true positives
points( true_positive$EWS_date
        , dat_df$cases[ which( dat_df$date %in% true_positive$EWS_date ) ]
        , col = "blue"
        , cex = 2
        #, pch = 16
        , lwd=2
)
#' Plot false positives
points( false_positive$EWS_date
        , dat_df$cases[ which( dat_df$date %in% false_positive$EWS_date ) ]
        , col = "orange"
        , cex = 2
        , pch = 16
        , lwd=2
)
#' Plot earliest true positive
points( true_positive_earliest
        , dat_df$cases[ which( dat_df$date %in% true_positive_earliest ) ]
        , col = "green"
        , cex = 2
        , pch = 16
        , lwd=2
)
points( wave_start_dates[ 1:7 ]
        , dat_df$cases[ which( dat_df$date %in% wave_start_dates ) ]
        , col = "red"
        , cex = 2
        , pch = 16
        , lwd=2
)
legend("topright"
       , legend = c("UK Covid-19 hospitalisations" , "Wave inflection points" , "Earliest true positive EWS" , "True positive EWS" , "False positive EWS")
       , col = c("black","red","green","blue","orange")
       , lty = c(1,NA,NA,NA,NA)
       , lwd = c(1.5,2,2,2,2) 
       , pch = c(NA,16,16,1,16)
       , cex = c(1.25,1.25,1.25,1.25,1.25)
)

#' Read in leading indicator data
#' Need to add +7 days for assumption on hospitalisation reporting (only for hospitalisations)
#' Read in leading indicator data
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_11/analysis_new_pango_stat_2/large_cluster_adjust_lgr_threshold_090/p_val_filter_no/dataframes_statistics')
filename = "tfps_pango_lineage_max_lgr.csv"
ews_base = fread( filename ) ;  ews_base[ , 1 ] = NULL ; ews_base = data.frame( ews_base )

#' Calculate robust Z-score time series (code taken from C:\Users\kdrake\OneDrive - Imperial College London\Documents\TFPS\tfps runs\2022_12\analysis\EWS_R_scripts/EWS_calc_threshold_lgr_th_090_pval_filter_no_HPC_array.R)
ews = ews_base[ , c( 1 , 2 ) ] # Select date column (1) and variable column (var_type)
colnames( ews ) <- c( "date" , "cases" ) # Rename columns
ews = ews[ !is.na( ews$cases ) , ]# Remove days with NA values
ews$date <- as.Date( ews$date , format = "%d/%m/%Y") # Change format of date
ews$time <- lubridate::decimal_date( ews$date ) # Add column for decimal date
ews$wday <- lubridate::wday( ews$date ) # Add column for day of week
ews = ews[ Reduce( '&' , lapply( ews , is.finite ) ) , ] # Remove days with Inf values
#' Add columns for running z-score and robust z-score (using median and MAD) calculating on add-one-in basis
running_z_score = data.frame()
running_z_score_robust = data.frame()
for ( n in 1 : dim( ews )[ 1 ] ){
  dat = ews[ 1 : n , 2 ]
  running_z_score[ n , 1 ] = ( ews[ n , 2 ] - mean( dat ) ) / sd( dat )
  #' Note - the distribution of variance of LGR is not always normal, 
  #' and there are some very high values for var(LGR) (e.g. outliers) which swamp the
  #' standard z-score, so it is better to use the robust z-score
  running_z_score_robust[ n , 1 ] = ( ews[ n , 2 ] - median( dat ) ) / mad( dat )
}
ews = cbind( ews , running_z_score, running_z_score_robust )
colnames( ews[ 5 ] ) <- c( "running_z_score" )
colnames( ews[ 6 ] ) <- c( "running_z_score_robust" )
#' plot robust Z-score time series
plot( ews[,1]
      , ews[ , 6 ] 
      , typ = "l" 
      , lwd = 1.5
      , col="black"
      , ylab="'Robust' Z-score of leading indicator"
      , xlab="Date"
      , cex.lab = 1.25
      , cex.axis = 1.25
      , xlim=c(min(dat_df$date),max(dat_df$date))
      , ylim = c(-5,15)
      #, ylim = c(min(na.omit(ews[,6])),max(na.omit(ews[,6])))
      , xaxs="i"
      , yaxs="i" #' Makes plot box tight to values
)#, ylim = c( 0 , 5 ) )
abline(h=0)
#' plot EWS threshold
abline(h=0.45,lty=2)

#' Plot true positives
points( true_positive$EWS_date
        , ews[ , 6 ][ which( ews[ , 1 ] %in% true_positive$EWS_date ) ]
        , col = "blue"
        , cex = 2
        #, pch = 16
        , lwd=2
)
#' Plot false positives
points( false_positive$EWS_date
        , ews[ , 6 ][ which( ews[ , 1 ] %in% false_positive$EWS_date ) ]
        , col = "orange"
        , cex = 2
        , pch = 16
        , lwd=2
)
#' Plot earliest true positive
points( true_positive_earliest
        , ews[ , 6 ][ which( ews[ , 1 ] %in% true_positive_earliest ) ]
        , col = "green"
        , cex = 2
        , pch = 16
        , lwd=2
)
#' Plot wave inflection/start points
points( wave_start_dates[ 1:7 ]
        , ews[ , 6 ][ which( ews[ , 1 ] %in% wave_start_dates ) ]
        , col = "red"
        , cex = 2
        , pch = 16
        , lwd=2
)

#' Plot leading indicator
plot( ews_base[ , 1 ] 
      , ews_base[ , 2 ] 
      , typ = "l"
      , lwd = 1.5
      , col="black"
      , ylab="Leading indicator"
      , xlab="Date"
      , xlim=c(min(dat_df$date),max(dat_df$date))
      , ylim = c(0,5)
      , cex.lab = 1.25
      , cex.axis = 1.25
      , xaxs="i"
      , yaxs="i" #' Makes plot box tight to values
)#, ylim = c( 0 , 5 ) )
#' Plot true positives
points( true_positive$EWS_date
        , ews_base[ , 2 ][ which( ews_base[ , 1 ] %in% true_positive$EWS_date ) ]
        , col = "blue"
        , cex = 2
        #, pch = 16
        , lwd=2
)
#' Plot false positives
points( false_positive$EWS_date
        , ews_base[ , 2 ][ which( ews_base[ , 1 ] %in% false_positive$EWS_date ) ]
        , col = "orange"
        , cex = 2
        , pch = 16
        , lwd=2
)
#' Plot earliest true positive
points( true_positive_earliest
        , ews_base[ , 2 ][ which( ews_base[ , 1 ] %in% true_positive_earliest ) ]
        , col = "green"
        , cex = 2
        , pch = 16
        , lwd=2
)
#' Plot wave inflection/start points
points( wave_start_dates[ 1:7 ]
        , ews_base[ , 2 ][ which( ews_base[ , 1 ]  %in% wave_start_dates ) ]
        , col = "red"
        , cex = 2
        , pch = 16
        , lwd=2
)


#############
#' ROC by wave and also by entire dataset[?]

#' Read in dataframe
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_11/analysis/EWS_df_combined")
wave_results_df = readRDS( "wave_results_df.rds" )

#' Initialise dataframe to record results of true and false positive EWS
ROC_df = data.frame(    "leading_indicator_type" = NA
                        , "TFPS_cluster_min_age" = NA
                        , "TFPS_cluster_max_age" = NA
                        , "TFPS_cluster_min_descendants" = NA
                        , "LGR_p_value_threshold" = NA
                        , "EWS_threshold" = NA
                        , "ews_tpr" = NA
                        , "ews_fpr" = NA
)

for ( lead_ind in filenames_prefix ){
  leading_indicator_type = substr( lead_ind , 1 , nchar( lead_ind ) - 1 )
  message( lead_ind )
  for ( min_age in c("07","14","28") ){
    for ( max_age in c("56","84") ){
      for ( min_desc in c(20,50,100,"proportional") ){
        for ( p_val_th in p_val_thresholds ){
          #' Loop through EWS threshold levels, which represent the classifier boundaries
          for ( ews_th in seq(0.25,5.00,0.25) ){
            
            #for ( i  in 1 : length( ews_thresholds ) ){
            #' Filter by leading indicator, scan variables, cluster filter (LGR p-value) and EWS threshold
            wave_results_df_filtered = subset( wave_results_df , wave_results_df$leading_indicator_type == leading_indicator_type &
                                                 wave_results_df$TFPS_cluster_min_age == min_age &
                                                 wave_results_df$TFPS_cluster_max_age == max_age &
                                                 wave_results_df$TFPS_cluster_min_descendants == min_desc &
                                                 wave_results_df$LGR_p_value_threshold == p_val_th &
                                                 wave_results_df$EWS_threshold == ews_th
            )
            #' Calculate True Positive (TPR) and False Positive Rates (FPR)
            ews_total = nrow( wave_results_df_filtered )
            ews_true_positive = nrow( subset( wave_results_df_filtered , wave_results_df_filtered$positive_EWS == TRUE ) )
            ews_false_positive = nrow( subset( wave_results_df_filtered , wave_results_df_filtered$positive_EWS == FALSE ) )
            ews_tpr = ews_true_positive / ews_total 
            #'**The all negatives value may need adapting given that need X consecutive values above EWS threshold positives before signal recorded i.e. do we include all the X data points above threshold that contribute to a positive EWS in the positive or negative results?**
            ews_fpr = ews_false_positive / ( length( ews$date ) - ews_total ) #' Number of incorrect positives out of all negatives. 
            #'** Could also calculate TPR and FPR per wave and calculate stats on lead/lag times**
            #'** Would also be good to record the earliest EWS per wave - so can show that higher threshold reduces FPR but increases lag time**
            
            #' And compile in dataframe
            new_row = ROC_df[ 1 , ]
            new_row$leading_indicator_type = leading_indicator_type
            new_row$TFPS_cluster_min_age = min_age
            new_row$TFPS_cluster_max_age = max_age
            new_row$TFPS_cluster_min_descendants = min_desc
            new_row$LGR_p_value_threshold = p_val_th
            new_row$EWS_threshold = ews_th
            new_row$ews_tpr = ews_tpr
            new_row$ews_fpr = ews_fpr
            
            ROC_df = rbind( ROC_df , new_row )
          }
        }
      }
    }
  }
}

View( ROC_df )

#' Plot ROC results for var(LGR) leading indicator (sample variance)
ROC_df_plot = subset( ROC_df , leading_indicator_type == "tfps_vlgr_samp" &
                        TFPS_cluster_min_age == "14" &
                        TFPS_cluster_max_age == "84" &
                        TFPS_cluster_min_descendants == "100" &
                        LGR_p_value_threshold == 10000
)

plot(ROC_df_plot[,8],ROC_df_plot[,7],typ="l",ylab="True Positive Rate (TPR)", xlab="False Positive Rate (FPR)")


library(ROCR)
data(ROCR.simple)
df <- data.frame(ROC_df_plot[,c(7,8)])
pred <- prediction(df$predictions, df$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)


#' Calculate Youden Index


#' Calculate Area Under Curve