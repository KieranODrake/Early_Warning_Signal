#' Generate EWS from TFPS data using same method applied to non_TFPS derived leading indicators
#' (this is different from the main analysis method used for the TFPS leading indicator data)
#' Purpose: Generate early warning signals (EWS) using leading indicator data 
#' A range of leading indicator thresholds are used for generating EWS. 
#' Absolute value of the leading indicator data is used.
#' Calculates various binary classifier stats for each leading indicator type.
#' See 'EWS_gen_time_window_TFPS_HPC.R' for calculation of ROC AUC only
#' But note Chicco & Jurman BioData Mining (2023) 16:4  https://doi.org/10.1186/s13040-023-00322-4
#' Author: Kieran Drake
#' Date: May 2023

#' 1 - Load data for each leading indicator type and format in dataframes with date and different variations in columns
#' 2 - Compute 'robust' z-scores for all leading indicator data
#' 3 - Compute whether z-score above or below EWS threshold
#' 4 - Load waves definitions: start/inflection date, Rt critical transition dates, wave band dates
#' 5 - For each leading indicator, wave, and EWS threshold, want to know number of TP, FP, TN and FN as well as earliest TP.
#' 6 - Can then calculate AUC and various other ROC stats (https://en.wikipedia.org/wiki/Receiver_operating_characteristic)
#'     for each wave and for whole period (average or total TP, FP, TN, FN across all waves).

library( magrittr )
library( lubridate )
#library( ggplot2 )
library( stringr )
library( zoo )
library( data.table )
library( sqldf )

#' Take array input on High Performance Computing (HPC)
params_run <- commandArgs( trailingOnly = TRUE )
params <- strsplit( params_run , "," ) #' multiple parameters input using array
#parent_sub_lgr_threshold = as.numeric( params[[1]][1] ) #' 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 false
#p_val_thresholds = as.numeric( params[[1]][2] ) #10000 #, 0.05 , 0.01 
parent_sub_lgr_threshold_char = as.character( params[[1]][1] ) #' 060 065 070 075 080 085 090 095 false
p_val_thresholds_char = as.character( params[[1]][2] ) #no # 005 , 001 

pslt = parent_sub_lgr_threshold_char
pvt = p_val_thresholds_char
###' 1 - Load data
#setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_060/p_val_filter_001/dataframes_statistics/')
filenames <- c(   "tfps_vlgr_samp.csv"      , "tfps_vlgr_simple_samp.csv"    , "tfps_vlgr_gam_samp.csv" 
                  , "tfps_vlgr_pop.csv"     , "tfps_vlgr_simple_pop.csv"     , "tfps_vlgr_gam_pop.csv"
                  , "tfps_vlgr_wtd.csv"
                  , "tfps_lgr_max.csv"      , "tfps_lgr_simple_max.csv"      , "tfps_lgr_gam_max.csv"
                  , "tfps_lgr_mean.csv"     , "tfps_lgr_simple_mean.csv"     , "tfps_lgr_gam_mean.csv"
                  , "tfps_lgr_wtd_mean.csv" , "tfps_lgr_simple_wtd_mean.csv" , "tfps_lgr_gam_wtd_mean.csv"
                  , "tfps_clock_outlier_max.csv" , "tfps_clock_outlier_mean.csv"
                  , "tfps_pango_lineage_max_lgr.csv"
)
#filenames_prefix <- paste( data.frame( str_split( filenames , pattern = ".csv" ) )[ 1 , ] , "_" , sep = "" )
filenames_prefix <- paste( data.frame( str_split( filenames , pattern = ".csv" ) )[ 1 , ] , sep = "" )
combined_df = data.frame()
for(i in 1:length(filenames)){
  temp <- fread( filenames[i] )
  temp[,1] = NULL
  temp = data.frame( temp )
  #new_name = paste0( filenames_prefix[ i ],"_df" )
  for(j in 2 : ncol( temp ) ){
    #colnames( temp )[ j ] <- paste0( colnames( temp )[ j ] , "_" , new_name )
    colnames( temp )[ j ] <- paste0( colnames( temp )[ j ] , "_" , filenames_prefix[ i ] )
  }
  #assign( new_name , temp )
  if( i == 1 ){
    combined_df <- temp  
  } else if ( i > 1 ){
    combined_df <- cbind( combined_df , temp[ , 2 : ncol(temp) ] )
  }
}
rm(temp)

#' Because the phylo trees were not published daily the TFPS time series have missing dates.
#' Therefore need to interpolate to generate dataset of daily datapoints.
#' Create list of dates in between weekly measurements
tfps_min_date = min( combined_df[,1] ) ; tfps_max_date = max( combined_df[,1] )
tfps_cal_dates = seq( tfps_min_date , tfps_max_date , 1 )
'%ni%' = Negate('%in%')
tfps_cal_dates_missing = tfps_cal_dates[ ( tfps_cal_dates %ni% combined_df[,1] ) ]
tfps_cal_dates_missing_df = data.frame(  "date" = tfps_cal_dates_missing )
#' Add columns for TFPS values
for (i in 2 : ncol( combined_df ) ){
  tfps_cal_dates_missing_df[ , i ] <- NA
  names( tfps_cal_dates_missing_df )[ i ] <- names( combined_df )[ i ]
}
#' Merge rows with actual data and empty rows for dates in between
combined_interpolate_df = merge( combined_df , tfps_cal_dates_missing_df 
                              #, by = intersect( names( comix_df ) , names( comix_cal_dates_missing_df ) )
                              , by.x = names( combined_df ) , by.y = names( tfps_cal_dates_missing_df )
                              , all.x = TRUE, all.y = TRUE)
#' Interpolate values in between actual data
combined_na_approx = combined_interpolate_df
for (i  in 2 : ncol( combined_interpolate_df ) ){
  #message(i)
  row_first = which( !is.na( combined_interpolate_df[,i] ) )[ 1 ] #' Find which row is the first to have a value
  row_last = max( which( !is.na( combined_interpolate_df[,i] ) ) )
  combined_na_approx[ row_first:row_last , i ] = na.approx( combined_interpolate_df[ row_first:row_last , i ] )
}
#' Check interpolation
#plot( combined_na_approx[,1],combined_na_approx[,2])
#points( combined_df[,1],combined_df[,2],col="red",pch=16)
#' Remove intermediate dataframes
rm(tfps_cal_dates_missing_df,combined_interpolate_df,tfps_max_date,tfps_min_date,tfps_cal_dates,tfps_cal_dates_missing)


##' Trim time series to maximum time period which all leading indicators cover
#' (the shortest time period happens to be one of the non-genomic leading indicators)
date_min = as.Date("2020-08-14")
date_max = as.Date("2022-03-09")
trim_combined_na_approx = subset( combined_na_approx , ( ( combined_na_approx$date >= date_min ) & ( combined_na_approx$date <= date_max ) ) )

#' Replace Inf values with NA and interpolate to replace all NAs 
tfps_df <- combined_df ; rm(combined_df)
for(i in 2 : ncol( tfps_df ) ){
  inf_index = which( is.infinite( tfps_df[ , i ] ) )
  tfps_df[ inf_index , i ] <- NA
  tfps_df[ , i ] = na.approx( tfps_df[ , i ] )
}


##' 2 - Compute 'robust' z score for TFPS leading indicators
#' Loop through time series for each set of variables in the leading indicator file
tfps_z_df = tfps_df
tfps_z_df[ , 2 : ncol( tfps_df ) ] = NA
for ( var_type in 2 : ncol( tfps_df ) ){
  #' Process EWS time series
  ews = tfps_df[ , c( 1 , var_type ) ] # Select date column (1) and variable column (var_type)
  colnames( ews ) <- c( "date" , "cases" ) # Rename columns
  ews$date <- as.Date( ews$date , format = "%d/%m/%Y") # Change format of date
  ews$time <- lubridate::decimal_date( ews$date ) # Add column for decimal date
  ews$wday <- lubridate::wday( ews$date ) # Add column for day of week
  #' Add columns for running z-score and robust z-score (using median and MAD) calculating on add-one-in basis
  running_z_score_robust = data.frame()
  for ( n in 1 : dim( ews )[ 1 ] ){
    dat = ews[ 1 : n , 2 ]
    #' Note - the distribution of variance of LGR is not always normal, 
    #' and there are some very high values for var(LGR) (e.g. outliers) which swamp the
    #' standard z-score, so it is better to use the robust z-score
    running_z_score_robust[ n , 1 ] = ( ews[ n , 2 ] - median( dat ) ) / mad( dat ) #running_z_score[ n , 1 ] = ( ews[ n , 2 ] - mean( dat ) ) / sd( dat )
  }
  ews = cbind( ews , running_z_score_robust )
  names( ews )[ 5 ] <- c( "running_z_score_robust" )
  tfps_z_df[ , var_type ] <- ews$running_z_score_robust
  #message( paste0(var_type," ",names( tfps_df )[ var_type ] )
}   
#plot( tfps_z_df[ , 1 ] , rowMeans( tfps_z_df[ , 2:ncol( tfps_z_df ) ] ) , ylim=c(-5,+5) , typ="l")
#abline(h=0)
#plot(tfps_z_df[,1],tfps_z_df[,2],ylim=c(-50,+30),typ="l")
#for(i in 3:ncol( tfps_z_df )){
#  lines(tfps_z_df[,1],tfps_z_df[,i],col=rainbow(ncol( tfps_z_df ))[i])
  #invisible( readline( prompt="Press [enter] to continue" ) )
#}
#for(i in 2:ncol( tfps_z_df )){
#  plot(tfps_z_df[,1],tfps_z_df[,i],col=rainbow(ncol( tfps_z_df ))[i],ylim=c(-50,30),main=names(tfps_z_df)[i],typ="l",lwd=2)
#  lines(tfps_df[,1],tfps_df[,i],col=rainbow(ncol( tfps_z_df ))[i+5],ylim=c(-50,30),main=names(tfps_z_df)[i],typ="l")
#  abline(h=0)
  #invisible( readline( prompt="Press [enter] to continue" ) )
#}
#plot(tfps_z_df[,1],tfps_z_df[,2],ylim=c(-30,+30),typ="l")
#abline(h=0)

#' 3 - Compute whether z-score above or below EWS threshold
#' Create an array. X columns for different leading indicators, 563 rows for dates and 201 z-dim for EWS threshold levels (-5 to +5 in 0.05 increments)
tfps_z_ews_array = array(data=NA
                         ,dim=c( nrow( tfps_z_df ) #' dates
                                 ,ncol( tfps_z_df ) - 1 #' leading indicator 'robust' z-score values -ignore first column which is the dates
                                 ,length( seq(-5,+5,0.05) ) #' Whether above or below EWS threshold
                         )
)
ews_count = 0
for( ews_th in seq(-5,5,0.05) ){
  ews_count = ews_count + 1 
  temp = !( tfps_z_df[ , 2 : ncol( tfps_z_df ) ] <= ews_th ) #' TRUE if above EWS threshold and FALSE if below
  #' Replace Inf and NA values with FALSE 
  for(i in 1 : ncol( temp ) ){
    inf_na_index = which( is.infinite( temp[ , i ] ) | is.na( temp[ , i ] ) )
    temp[ inf_na_index , i ] <- FALSE
  }
  tfps_z_ews_array[ , , ews_count ] = temp
}

##' 4 - Load waves definitions: start/inflection date, Rt critical transition dates, wave band dates
#' Wave start dates (average where a wave has more than one start date as 
#' produced from the 12 models incorporating 'low resolution wave filter' - see
#' 21 June 2022 update slides)
wave_start_dates <- c(                        #' 1 Wuhan
                       as.Date("2020-08-19") #' 2 B.1.177
                      ,as.Date("2020-11-29") #' 3 Alpha
                      ,as.Date("2021-05-11") #' 4 Delta
                      ,as.Date("2021-08-03") #' 5 Delta
                      ,as.Date("2021-09-27") #' 6 Delta
                      ,as.Date("2021-11-26") #' 7 Omicron
                      ,as.Date("2022-02-21") #' 8 Omicron
                      #' 9 Omicron
)

Rt_critical_transitions <- c(                        #' 1 Wuhan
                               as.Date("2020-09-06") #' 2 B.1.177
                              ,as.Date("2020-12-13") #' 3 Alpha
                              ,as.Date("2021-05-22") #' 4 Delta (1)
                              ,as.Date("2021-08-19") #' 5 Delta (2)
                              ,as.Date("2021-10-13") #' 6 Delta (3)
                              ,as.Date("2021-11-25") #' 7 Omicron BA.1
                              ,as.Date("2022-03-11") #' 8 Omicron BA.2
                                                     #' 9 Omicron
)

#' Alternative method for defining wave bands as some EWS dates not falling within bands
#' Manually chosen date roughly half way between peak and trough on the wave decline  
wave_bands_start <- c(#as.Date("2020-01-01"), # 1  Wuhan
                      as.Date("2020-05-01"), # 2  B.1.177
                      as.Date("2020-11-23"), # 3  Alpha
                      as.Date("2021-02-01"), # 4  Delta
                      as.Date("2021-08-01"), # 5  Delta
                      as.Date("2021-09-18"), # 6  Delta
                      as.Date("2021-11-11"), # 7  Omicron
                      # as.Date("2022-01-05"), # 8  Omicron Not an obvious separate peak in hospitalisation data
                      as.Date("2022-01-29") # 9  Omicron
                      #as.Date("2022-04-23")  # 10 Omicron
)
wave_bands_end <- c(#as.Date("2020-04-30"), # 1  Wuhan
                      as.Date("2020-11-22"), # 2  B.1.177
                      as.Date("2021-01-31"), # 3  Alpha
                      as.Date("2021-07-30"), # 4  Delta
                      as.Date("2021-09-17"), # 5  Delta
                      as.Date("2021-11-10"), # 6  Delta
                      as.Date("2022-01-28"), # 7  Omicron
                      #as.Date("2022-01-18"), # 8  Omicron Not an obvious separate peak in hospitalisation data
                      as.Date("2022-04-22") # 9  Omicron
                      #as.Date("2022-07-30")  # 10 Omicron
)

#' Wave table
wave_df = data.frame(   "wave_n" = seq(2,8,1)
                        , "wave_name" = c("Unk","Alpha","Delta_1", "Delta_2","Delta_3","Omicron_1","Omicron_2")
                        , "pango" = c("B.1.177","B.1.1.7","B.1.617.2","B.1.617.2","B.1.617.2","BA.1","BA.2")
                        , "wave_inflection_start_date" = wave_start_dates
                        , "wave_band_start" = wave_bands_start 
                        , "wave_band_end" = wave_bands_end
) #possibly AY.4 instead of B.1.617.2

#' Define periods when value above threshold considered TRUE EWS
true_ews_start = -30 #' t-30 as per Proverbio et al - used to avoid interference from previous wave 
true_ews_end = + 10 #' t+10 to account for some later earliest True Positie EWS in the TFPS data. Compares with t+5 as per Proverbio et al 
true_start_dates = list() ; true_rt_dates = list()
for(i in 1:7){
  #' Dates where value above threshold would be a TRUE Positive (TP) EWS and below threshold would be a FALSE negative (FN) EWS
  true_start_dates[[i]] = seq( wave_start_dates[i] + true_ews_start , wave_start_dates[i] + true_ews_end , 1 )
  true_rt_dates[[i]] = seq( Rt_critical_transitions[i] + true_ews_start , Rt_critical_transitions[i] + true_ews_end , 1 )
}
#' On dates other than those defined above, the classification is reversed:
#' False positive (FP) if above EWS threshold and TRUE negative (TN) if below EWS threshold
#' Apply to array where TRUE simply indicates above the EWS threshold and FALSE indicates below or equal to EWS threshold
true_start_dates_index = which( tfps_z_df[,1] %in% unlist(true_start_dates) )
true_rt_dates_index = which( tfps_z_df[,1] %in% unlist(true_rt_dates) )

tfps_z_ews_sd_classification_array = tfps_z_ews_array #' Initialise new array for classifying signals using existing array as template
tfps_z_ews_rt_classification_array = tfps_z_ews_array #' Initialise new array for classifying signals using existing array as template

for( k in 1 : nrow( tfps_z_df ) ){ #' loop through dates
  for( m in 1 : ( ncol( tfps_z_df ) -1 ) ){ #' loop through leading indicators
    for( n in 1 : dim( tfps_z_ews_array )[ 3 ] ){ #' loop through EWS thresholds
      #' Classify based on wave start dates
      if( (k %in% true_start_dates_index ) & ( tfps_z_ews_array[k,m,n] == TRUE ) ){ #' WITHIN EWS wave range & ABOVE EWS threshold
        tfps_z_ews_sd_classification_array[k,m,n] = "TP"
      }
      if( (k %in% true_start_dates_index ) & ( tfps_z_ews_array[k,m,n] == FALSE ) ){ #' WITHIN EWS wave range & BELOW EWS threshold
        tfps_z_ews_sd_classification_array[k,m,n] = "FN"
      }
      if( (k %ni% true_start_dates_index ) & ( tfps_z_ews_array[k,m,n] == TRUE ) ){ #' OUTSIDE EWS wave range & ABOVE EWS threshold
        tfps_z_ews_sd_classification_array[k,m,n] = "FP"
      }
      if( (k %ni% true_start_dates_index ) & ( tfps_z_ews_array[k,m,n] == FALSE ) ){ #' OUTSIDE EWS wave range & BELOW EWS threshold
        tfps_z_ews_sd_classification_array[k,m,n] = "TN"
      }
      #' Classify based on Rt critical transition dates
      if( (k %in% true_rt_dates_index ) & ( tfps_z_ews_array[k,m,n] == TRUE ) ){ #' WITHIN EWS wave range & ABOVE EWS threshold
        tfps_z_ews_rt_classification_array[k,m,n] = "TP"
      }
      if( (k %in% true_rt_dates_index ) & ( tfps_z_ews_array[k,m,n] == FALSE ) ){ #' WITHIN EWS wave range & BELOW EWS threshold
        tfps_z_ews_rt_classification_array[k,m,n] = "FN"
      }
      if( (k %ni% true_rt_dates_index ) & ( tfps_z_ews_array[k,m,n] == TRUE ) ){ #' OUTSIDE EWS wave range & ABOVE EWS threshold
        tfps_z_ews_rt_classification_array[k,m,n] = "FP"
      }
      if( (k %ni% true_rt_dates_index ) & ( tfps_z_ews_array[k,m,n] == FALSE ) ){ #' OUTSIDE EWS wave range & BELOW EWS threshold
        tfps_z_ews_rt_classification_array[k,m,n] = "TN"
      }
    }
  }
}


#' 5 - For each leading indicator, wave, and EWS threshold, want to know number of TP, FP, TN and FN as well as earliest TP.
#'  classification types
#' 201 EWS thresholds
#' 7 waves + 1 total + 1 total ex.B.1.177

#' Create range of dates associated with each wave
date_index_list = list()
for(p in 1:7){
  date_index_temp = which( tfps_z_df[,1] %in% seq( wave_bands_start[ p ], wave_bands_end[ p ] , 1 ) )
  assign( paste0( "date_index_w" , p+1 ), date_index_temp )
  date_index_list[[p]] <- date_index_temp
  names( date_index_list )[p] <- paste0( "date_index_w" , p+1 )
}
#' Create date range for total ex. B.1.177
date_index_total <- which( tfps_z_df[,1] %in% seq( wave_bands_start[ 1 ], wave_bands_end[ 7 ] , 1 ) )
#' Create date range for total ex. B.1.177
date_index_total_ex_B.1.177 <- which( tfps_z_df[,1] %in% seq( wave_bands_start[ 2 ], wave_bands_end[ 7 ] , 1 ) )

#' Create date range for total
date_index_total            <- which( tfps_z_df[,1] %in% seq( wave_bands_start[ 1 ], wave_bands_end[ 7 ] , 1 ) )
#' Create date range for total ex. B.1.177 (not enough data)
date_index_total_ex_B.1.177 <- which( tfps_z_df[,1] %in% seq( wave_bands_start[ 2 ], wave_bands_end[ 7 ] , 1 ) )
#' Create date range for total ex. B.1.177 and BA.2 (not enough data)
date_index_total_ex_B.1.177_BA.1 <- which( tfps_z_df[,1] %in% seq( wave_bands_start[ 2 ], wave_bands_end[ 6 ] , 1 ) )
#' Add to list of lists
date_index_list[[ 8]] = date_index_total                 ; names(date_index_list)[ 8] <- "date_index_total"
date_index_list[[ 9]] = date_index_total_ex_B.1.177      ; names(date_index_list)[ 9] <- "date_index_total_ex_B.1.177"
date_index_list[[10]] = date_index_total_ex_B.1.177_BA.1 ; names(date_index_list)[10] <- "date_index_total_ex_B.1.177_BA.1"

#' **Previous section corresponds to '1_EWS_gen_time_window_non_TFPS_processing.R' for non-TFPS leading indicators**
#' **and following section corresponds to '2_EWS_gen_time_window_non_TFPS_multi_stat.R'**

#' 6 Calculation
#' 6.1 - Loop through the two different ways of calculating the time window (anchoring to Rt critical transition date or hospitalisation wave inflection date)
#' 6.2 - Loop through date ranges (individual waves, inc all waves, all waves exc. B.1.177 and all waves exc B.1.177 & BA.2)
#' 6.3 - Loop through leading indicators
#' 6.4 - Calculate ROC stats and AUC for each wave etc.

#' Define area under curve calculation (AUC) function using trapezoidal rule of integration
#' https://en.wikipedia.org/wiki/Trapezoidal_rule and https://stackoverflow.com/questions/7358738/how-to-calculate-integral-with-r
AUC <- function(x, y){
  sum(diff(x)*rollmean(y,2))
}

#ews_class_total = table( li_z_ews_rt_classification_array[ c( date_index_total ) , 1 , 201 ] )

for( arr in 1:2 ){
  message(arr)
  #' initialise output data.frame template
  ROC_stat_types = c("TP","FP","TN","FN","TPR","FPR","TNR","PPV","NPV","F1","FMI","MCC","normMCC","ROC_AUC")
  output = data.frame( matrix( NA , nrow = ( ( ncol( tfps_z_df ) - 1 ) * dim( tfps_z_ews_rt_classification_array )[ 3 ]) #' <number of leading indicators> * <number of ews threshold values>
#**MIGHT NEED MORE COLUMNS HERE FOR TFPS VS non-TFPS (to include min and max age and min descendants, pslt, pvt)**
                                  , ncol = 2 + ( length( ROC_stat_types ) * length( date_index_list ) ) #' number of date ranges plus 2 columns for the leading indicator type and the ews threshold value
  ) 
  ) 
  #' Add column names to output template
  colnames(output)[1:2] = c("leading_indicator_type","ews_threshold")
  counter=0
  for(wave_no in 2:(length(date_index_list)+1)){
    print( wave_no )
    for( ROC_stat in ROC_stat_types ){
      counter = counter + 1
      if       (wave_no<= 8){ prefix <- c(paste0("w",wave_no,"_")) 
      } else if(wave_no== 9){ prefix <- c(paste0("w_all_"       ))
      } else if(wave_no==10){ prefix <- c(paste0("w_all_ex_2_"  )) 
      } else if(wave_no==11){ prefix <- c(paste0("w_all_ex_2_8_")) }
      new_col_name = paste0( prefix , ROC_stat )
      message( new_col_name )
      colnames( output )[ counter+2 ] = new_col_name
    }
  }
  
  date_range_counter = 0
  for( date_index in date_index_list ){
    date_range_counter = date_range_counter + 1
    
    for( li in 1 : dim( tfps_z_ews_rt_classification_array )[ 2 ] ){
      #' Select between Rt critical transition date and wave start date
      if       (arr == 1){ temp_array = tfps_z_ews_rt_classification_array[ date_index , li ,  ] #' Using Rt critical transition as anchor point for time window
      } else if(arr == 2){ temp_array = tfps_z_ews_sd_classification_array[ date_index , li ,  ] }#' Using hospitalisation wave start date as anchor point for time window
      message( "Date range: ", date_range_counter , ". Leading indicator: ", li )
      
      #' Analyse EWS across all EWS thresholds for the particular date range (date_index) and leading indicator (li)
      temp_val_list = list()
      for(i in 1:dim( temp_array )[ 2 ]){
        temp_val = list( temp_array[  , i ] )
        temp_val_list = append( temp_val_list , temp_val )
      }
      ews_class_total = lapply( temp_val_list , table ) #( li_z_ews_rt_classification_array[ date_index , li , ews_th ] )
      out = data.frame()
      temp_output <- lapply( seq_along( ews_class_total ) 
                             , function( x , i ){   TP = as.numeric( ifelse( is.na( x[[i]]["TP"] ) , 0 , x[[i]]["TP"] ) ) #' TP
                             FP = as.numeric( ifelse( is.na( x[[i]]["FP"] ) , 0 , x[[i]]["FP"] ) ) #' FP
                             TN = as.numeric( ifelse( is.na( x[[i]]["TN"] ) , 0 , x[[i]]["TN"] ) ) #' TN
                             FN = as.numeric( ifelse( is.na( x[[i]]["FN"] ) , 0 , x[[i]]["FN"] ) ) #' FN
                             TPR = TP / (TP + FN) #' TPR True Positive Rate aka Sensitivity, true positive rate, recall, hit rate
                             FPR = FP / (FP + TN) #' FPR False Positive Rate = 1 - Specificity
                             TNR = TN / (TN + FP) #' TNR = True Negative Rate aka Specificity, selectivity, true negative rate
                             PPV = TP / (TP + FP)  #' PPV aka Precision, positive predictive value
                             NPV = TN / (TN + FN) #' NPV = Negative predictive value
                             F1 = (2*TP) / (2*TP + FP + FN) #' F1
                             FMI = sqrt( PPV * TPR )  #' Fowlkes-Mallows index
                             #' Need to break down MCC calculation as otherwise results in NA with message 'NAs produced by integer overflow'
                             #a = as.numeric(TP+FP) ; b = as.numeric(TP+FN) ; c = as.numeric(TN+FP) ; d = as.numeric(TN+FN)
                             #a = as.integer(TP+FP) ; b = as.integer(TP+FN) ; c = as.integer(TN+FP) ; d = as.integer(TN+FN)
                             #e = a*b ; f = e*c ; g = f*d
                             MCC = ((TP*TN) - (FP*FN)) / sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) #' combines the 4 basic rates of confusion matrix
                             #MCC = ((TP*TN) - (FP*FN)) / sqrt( a*b ) #' combines the 4 basic rates of confusion matrix
                             normMCC = (MCC + 1) / 2
                             #' Add ROC stats to dataframe
                             out[ 1, c(seq(1,13,1))] = c(TP,FP,TN,FN,TPR,FPR,TNR,PPV,NPV,F1,FMI,MCC,normMCC)
                             return(out)
                             } 
                             , x = ews_class_total )
      temp_output_df = data.frame()
      for( j in 1: length( temp_output ) ){
        temp_output_df = rbind( temp_output_df , as.data.frame( temp_output[[j]] ) )
      }
      
      #' Calculate area under curve (AUC) for TPR vs FPR
      #' Rows to include must be for a single leading indicator but across all EWS thresholds
      TPR_FPR = data.frame( "FPR" = temp_output_df[ ,  6 ] , "TPR" = temp_output_df[ ,  5 ] )
      names( TPR_FPR )[c(1,2)] <- c("FPR","TPR") ; TPR_FPR = TPR_FPR[ with( TPR_FPR , order( FPR , TPR ) ) , ]
      ROC_AUC = AUC( x = TPR_FPR$FPR, y = TPR_FPR$TPR )
      temp_output_df[ , length( ROC_stat_types ) ] = ROC_AUC
      colnames( temp_output_df ) = ROC_stat_types
      #' Add data for this date range and leading indicator to the main dataframe
      #' First need to determine location in the dataframe
      row_start = 1 + ( (li-1) * dim( tfps_z_ews_rt_classification_array )[ 3 ] )
      row_end   = dim( tfps_z_ews_rt_classification_array )[ 3 ] + ( (li-1) * dim( tfps_z_ews_rt_classification_array )[ 3 ] )
      col_start = 3 + ( ( date_range_counter -1 ) * length(ROC_stat_types) )
      col_end   = 16 + ( ( date_range_counter -1 ) * length(ROC_stat_types) )
      #' Check source and destination dimensions match
      dim( output[ row_start:row_end , col_start:col_end ] ) == dim( temp_output_df )
      #' Then add data to final dataframe
      output[ row_start:row_end , col_start:col_end ] = temp_output_df
      
      if(date_range_counter == 2){
        output[ row_start:row_end ,  1 ] = colnames(tfps_z_df)[li+1]
        output[ row_start:row_end ,  2 ] = seq(-5,5,0.05)#[ews_th]
      }
    } #' li for loop end
  } #' date_index for loop end
  
  #' Change name of output dataframe depending on which dataset was used
  if       (arr == 1){ new_name = "output_rt" #; assign( output , output_rt ) } #' Using Rt critical transition as anchor point for time window
  } else if(arr == 2){ new_name = "output_sd" }#; assign( output , output_sd ) }#' Using hospitalisation wave start date as anchor point for time window
  assign( new_name , output )
  
}

#' Save files
saveRDS( output_sd , file = paste0( "ROC_multi_stat_sd_df_lgr_th_",pslt,"_pval_th_",pvt,"_10d.rds" ) )
saveRDS( output_rt , file = paste0( "ROC_multi_stat_rt_df_lgr_th_",pslt,"_pval_th_",pvt,"_10d.rds" ) )
