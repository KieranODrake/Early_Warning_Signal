#' Splitting the counting of false positive EWS between those generated before 
#' the earliest true positive and those generated after the earliest true 
#' positive for each parameter set
#' 

library(data.table)

#' Load results dataframes
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined")
wave_results_df = readRDS( "wave_results_df.rds" )
wave_results_analysis_reshaped = readRDS( "wave_results_analysis_reshaped.rds" )

message("wave_results_df loaded")
#' Add column to wave_results_df with single identifier for leading indicator parameter set 
wave_results_df$variables_code = paste0(  wave_results_df$leading_indicator_type, "_"
                                          , wave_results_df$TFPS_cluster_min_age , "_"
                                          , wave_results_df$TFPS_cluster_max_age , "_"
                                          , wave_results_df$TFPS_cluster_min_descendants , "_"
                                          , wave_results_df$LGR_p_value_threshold , "_"
                                          , wave_results_df$parent_sub_LGR_threshold , "_"
                                          , wave_results_df$EWS_threshold
)

#' Link ews_results_df to earliest true positive dates in wave_results_analysis_df
earliest_tp_df = data.frame(   "variables_code" = wave_results_analysis_reshaped$variables_code
                               , wave_results_analysis_reshaped[ , c( seq(9,27,3) ) ]
                               , "w2_n_tp" = replicate(nrow(wave_results_analysis_reshaped),NA), "w2_n_fp_before_tp" = replicate(nrow(wave_results_analysis_reshaped),NA) , "w2_n_fp_after_tp" = replicate(nrow(wave_results_analysis_reshaped),NA)
                               , "w3_n_tp" = replicate(nrow(wave_results_analysis_reshaped),NA), "w3_n_fp_before_tp" = replicate(nrow(wave_results_analysis_reshaped),NA) , "w3_n_fp_after_tp" = replicate(nrow(wave_results_analysis_reshaped),NA)
                               , "w4_n_tp" = replicate(nrow(wave_results_analysis_reshaped),NA), "w4_n_fp_before_tp" = replicate(nrow(wave_results_analysis_reshaped),NA) , "w4_n_fp_after_tp" = replicate(nrow(wave_results_analysis_reshaped),NA)
                               , "w5_n_tp" = replicate(nrow(wave_results_analysis_reshaped),NA), "w5_n_fp_before_tp" = replicate(nrow(wave_results_analysis_reshaped),NA) , "w5_n_fp_after_tp" = replicate(nrow(wave_results_analysis_reshaped),NA)
                               , "w6_n_tp" = replicate(nrow(wave_results_analysis_reshaped),NA), "w6_n_fp_before_tp" = replicate(nrow(wave_results_analysis_reshaped),NA) , "w6_n_fp_after_tp" = replicate(nrow(wave_results_analysis_reshaped),NA)
                               , "w7_n_tp" = replicate(nrow(wave_results_analysis_reshaped),NA), "w7_n_fp_before_tp" = replicate(nrow(wave_results_analysis_reshaped),NA) , "w7_n_fp_after_tp" = replicate(nrow(wave_results_analysis_reshaped),NA)
                               , "w8_n_tp" = replicate(nrow(wave_results_analysis_reshaped),NA), "w8_n_fp_before_tp" = replicate(nrow(wave_results_analysis_reshaped),NA) , "w8_n_fp_after_tp" = replicate(nrow(wave_results_analysis_reshaped),NA)
)

#' Prepare list for lapply. Separate lists for each wave
earliest_tp_dt = as.data.table( earliest_tp_df )
rm(earliest_tp_df)
#' Reduce wave_results_df to essential columns for subsetting and optimise ability to subset 
#' https://stackoverflow.com/questions/14139586/how-to-optimize-subsetting-from-a-large-dataset
wave_results_dt = as.data.table( wave_results_df[ , c(74,66,73,9) ] )
wave_results_dt = subset( wave_results_dt , !is.na( wave_results_dt$wave_n ) )
wave_results_dt = wave_results_dt[ order( variables_code ) ]
wave_results_dt = setkey( wave_results_dt , variables_code )

#' remove dataframes to save memory now that smaller data tables have been created
rm( wave_results_df , wave_results_analysis_reshaped )
gc()

#' Create function to split the false positive EWS into before and after the earliest true positive
#' @param wave_results_dt Data table of wave results created using '1_false_positive_split_prep_for_HPC.R' and saved in 'false_positive_split_env_HPC_start.RData'.
#' @param earliest_tp_dt Predefined data table to be populated
#' @param wave_number the number of the wave so that filters can be applied and columns selected
#' @return Dataframe populated with false positives split between those before and those after the earliest true positive

fp_split <- function( wave_results , earliest_tp , wave_number ) 
{
  #' Split wave results into true and false positives
  wave_results_false_dt = subset( wave_results , wave_results$wave_n == wave_number & wave_results$positive_EWS == FALSE )
  wave_results_true_dt = subset( wave_results , wave_results$wave_n == wave_number & wave_results$positive_EWS == TRUE )
  rm( wave_results )
  gc()   
  
  #' Select earliest true positive list appropriate to wave
  if       ( wave_number == 2 ){ earliest_tp = earliest_tp[, c(1,2, 9,10,11) ]
  } else if( wave_number == 3 ){ earliest_tp = earliest_tp[, c(1,3,12,13,14) ]
  } else if( wave_number == 4 ){ earliest_tp = earliest_tp[, c(1,4,15,16,17) ]
  } else if( wave_number == 5 ){ earliest_tp = earliest_tp[, c(1,5,18,19,20) ]
  } else if( wave_number == 6 ){ earliest_tp = earliest_tp[, c(1,6,21,22,23) ]
  } else if( wave_number == 7 ){ earliest_tp = earliest_tp[, c(1,7,24,25,26) ]
  } else if( wave_number == 8 ){ earliest_tp = earliest_tp[, c(1,8,27,28,29) ]}
  #' Split data table into list of lists
  earliest_tp_list <- split( earliest_tp , seq( nrow( earliest_tp ) ) ) #' splits by row
  earliest_tp_list <- lapply( earliest_tp_list , as.list ) #' Turns into list of lists instead of list of single row dataframes
  rm( earliest_tp )
  gc()
  
  #' Run analysis
  output_list <- lapply( seq_along( earliest_tp_list ) , function( x , wf, wt, i){
    id = x[[i]]$variables_code
    #' Filter true and false positives by individual parameter set (id)
    wf = setkey( wf , variables_code )    
    temp_wf = wf[ J( id ) ]
    wt = setkey( wt , variables_code )    
    temp_wt = wt[ J( id ) ]
    
    #' Count the number of false positives before the earliest true positive, false positives after earliest true positive and number of true positives
    if ( (sum( !is.na( temp_wf$wave_n ) ) > 0) &  (nrow(temp_wf) > 0) ){ #' Set to NA if only NAs in the wave column
      x[[i]][4] = nrow( subset( temp_wf , temp_wf$EWS_date <  x[[i]][2] & !is.na(temp_wf$wave_n ) ) )
      x[[i]][5] = nrow( subset( temp_wf , temp_wf$EWS_date >= x[[i]][2] & !is.na(temp_wf$wave_n ) ) )
    } else {   
      x[[i]][4] = NA
      x[[i]][5] = NA
    }
    if ( ( sum( !is.na( temp_wt$wave_n ) ) > 0) &  ( nrow(temp_wt) > 0 ) ){ #' Set to NA if only NAs in the wave column
      x[[i]][3] = nrow( subset( temp_wt ,  !is.na(temp_wt$wave_n ) ) )
    } else {   
      x[[i]][3] = NA
    }
    
    return( x[[i]] ) 
    
  } , x = earliest_tp_list , wf = wave_results_false_dt , wt = wave_results_true_dt )
  
  
  return( output_list ) # Return from function
}

#' Save environment for use in High Performance Computing (HPC) cluster
save.image( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/false_positive_split_2023_03/false_positive_split_env_HPC_start_v2.RData" )
