#' Assessing sensitivity of SARS-CoV-2 early warning signal analysis, using 
#' the Transmission Fitness Polymorphism (TFP) Scanner, to region input into 
#' TFP Scannner (adm1 vs adm2) and generation time (Tg) (6.5 days vs 5 days)
#' Author: Kieran Drake
#' Date: October 2023
#' 

#' Summary
#' 1: Compare number of sequences used with adm1 and adm2 analyses (some sequences don't have adm2 value)
#' 1a: Generate time series of number of sequences in adm1 and adm2 data sets
#' 1b: Compute percentage difference
#' 1c: Plot time series (from 1a and 1b) along with UK SC2 hospitalisations
#' 2: Compare leading indicator time series created under different region type (adm1 vs adm2) analyses
#' 3: Compare leading indicator time series created with different generation times (6.5 days and 5 days)
#' 4: Extract data on lead times and false positives


#### 1 Sample frequency ########################################################
#' 1a

#' Read in sequence metadata associated with adm1 analysis (in preprint)
load("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/tfps-env-HPC-start_v5.RData" )
amd_adj_adm_1 = amd_adj
rm( amd_adj , bad_dates, no_sample_date , samples_to_be_removed_unique_name , scan_error_list , tre_list , tre_no_tips_in_meta_list, tre_read_error_list)

#' Read in sequence metadata associated with adm2 analysis
load("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/tfps-env-HPC-start_v5_adm2.RData" )
amd_adj_adm_2 = amd_adj
rm( amd_adj , bad_dates, no_sample_date , samples_to_be_removed_unique_name , scan_error_list , tre_list , tre_no_tips_in_meta_list, tre_read_error_list)

#' Create time series of number of sequences and % difference in number
library(data.table)
adm_1_seq_ts_df = as.data.frame( table(amd_adj_adm_1$sample_date) )
colnames( adm_1_seq_ts_df )[1] <- "sample_date"
adm_1_seq_ts_df$sample_date = as.Date( adm_1_seq_ts_df$sample_date )
adm_2_seq_ts_df = as.data.frame( table(amd_adj_adm_2$sample_date) )
colnames( adm_2_seq_ts_df )[1] <- "sample_date"
adm_2_seq_ts_df$sample_date = as.Date( adm_2_seq_ts_df$sample_date )
seq_n_ts <- merge( x = adm_1_seq_ts_df , y = adm_2_seq_ts_df
                               , by.x = "sample_date" , by.y = "sample_date"
                               , all.x = TRUE 
                               , all.y = FALSE )
colnames( seq_n_ts )[2] <- "sample_freq_adm1" ; colnames( seq_n_ts )[3] <- "sample_freq_adm2"
#' Replace NAs with 0
seq_n_ts$sample_freq_adm2[ is.na(seq_n_ts$sample_freq_adm2)] <- 0

#' 1b
#' Compute difference in sample frequency between adm1 and adm2 datasets
seq_n_ts$difference_abs = seq_n_ts$sample_freq_adm1 - seq_n_ts$sample_freq_adm2
seq_n_ts$difference_perc = (seq_n_ts$sample_freq_adm2 / seq_n_ts$sample_freq_adm1 -1)

#' 1c
#' Load case/hospitalisation data
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs'
#filename <- "data_2022-May-09 - UK Covid-19 cases by sample date.csv"
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
#dat_type <- "cases"
dat_type <- "hospitalisations"
#' data_load() @ C:\Users\kdrake\GitHub\Early_Warning_Signal\6c
library(data.table)
dat_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)
par( mfrow = c( 2 , 1 ) )
plot( dat_df$date
      , dat_df$cases
      , xlab = ""
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

#' Plot of sample frequency (adm1 and adm2) and hospitalisations - absolute values
plot( seq_n_ts$sample_date , seq_n_ts$sample_freq_adm1, lty=1 , xlab="", ylab="sample frequency",
      , cex.lab = 1.25 , cex.axis = 1.25 , lwd = 1.5
      , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) 
      ) 
lines( seq_n_ts$sample_date , seq_n_ts$sample_freq_adm1, col = "black")
lines( seq_n_ts$sample_date , seq_n_ts$sample_freq_adm2, col = "red")
legend("topleft",legend=c("adm1","adm2"),col=c("black","red"),lty=c(1,1))

#' Plot of sample frequency (adm1 and adm2) and hospitalisations - absolute difference
plot( seq_n_ts$sample_date , seq_n_ts$difference_abs, lty=1 , xlab="", ylab="sample frequency difference (adm1-adm2)",
      , cex.lab = 1.25 , cex.axis = 1.25 , lwd = 1.5
      , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) 
) 
lines( seq_n_ts$sample_date , seq_n_ts$difference_abs, col = "black")

#' Plot of sample frequency (adm1 and adm2) and hospitalisations - percentage difference
plot( seq_n_ts$sample_date , seq_n_ts$difference_perc*100, lty=1 , xlab="", ylab="sample frequency difference (% change)",
      , cex.lab = 1.25 , cex.axis = 1.25 , lwd = 1.5 , col="green"
      , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) 
) 
lines( seq_n_ts$sample_date , seq_n_ts$difference_perc *100 , col = "green")

#### 2 Leading indicator time series  ############################################
#' 2a: Read in leading indicator time series and calculate 'robust' z-score
#' Code below adapted from C:\Users\kdrake\OneDrive - Imperial College London\Documents\TFPS\tfps runs\2023_02\quantile_historical\analysis\HPC_scripts\EWS_calc_threshold_lgr_th_085_pval_filter_001_HPC_array.R
filename_1_adm2 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm2_Tg6_5/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_001/dataframes_statistics/tfps_pango_lineage_max_lgr.csv"
filename_1_adm1 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_001/dataframes_statistics/tfps_pango_lineage_max_lgr.csv"
filename_2_adm2 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm2_Tg6_5/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_no/dataframes_statistics/tfps_lgr_simple_mean.csv"
filename_2_adm1 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_no/dataframes_statistics/tfps_lgr_simple_mean.csv"

for ( filename in c( filename_1_adm1 , filename_1_adm2 , filename_2_adm1 , filename_2_adm2 ) ){

  ews_base = data.table::fread( filename ) 
  
  ews_base[ , 1 ] = NULL ; ews_base = data.frame( ews_base )
  
  #' Select relevant column for parameter set required
  if( filename %in% c( filename_1_adm2, filename_1_adm1 ) ){ 
    ews = ews_base[ , c("date","mina07_maxa56_md020")] #' Only need data for one set of parameters for this part
  } else if( filename %in% c(filename_2_adm2, filename_2_adm1) ){
    ews = ews_base[ , c("date","mina14_maxa56_md020")] 
  }
  
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
  colnames( ews )[ 5 ] <- c( "running_z_score" )
  colnames( ews )[ 6 ] <- c( "running_z_score_robust" )
  if( filename == filename_1_adm1 ) { filename_1_adm1_ews = ews } 
  if( filename == filename_1_adm2 ) { filename_1_adm2_ews = ews } 
  if( filename == filename_2_adm1 ) { filename_2_adm1_ews = ews } 
  if( filename == filename_2_adm2 ) { filename_2_adm2_ews = ews } 
}
  
  
  #' Four plots: hopsitalisations, leading indicator value, z-score, z-score difference
  par( mfrow = c( 4 , 1 ) )
  #' plot hospitalisations
  plot( x = dat_df$date, y = dat_df$cases
        , xlab = "" , ylab = "UK Covid-19 hospitalisations"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )  , ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  #' leading indicator value
  plot( x = filename_1_adm1_ews$date , y = filename_1_adm1_ews$cases
        , xlab = ""  , ylab = "leading indicator value"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_1_adm2_ews$date , y = filename_1_adm2_ews$cases
        , xlab = ""  , ylab = "leading indicator value"
        , cex.lab = 1.25  , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "red" #data_colour[2]
  )
  legend("topleft",legend=c("adm1","adm2"),col=c("black","red"),lty=c(1,1))
  
  #' robust z-score for leading indicator
  plot( x = filename_1_adm1_ews$date  , y = filename_1_adm1_ews$running_z_score_robust
        , xlab = ""  , ylab = "'Robust' z-score of leading indicator"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_1_adm2_ews$date , y = filename_1_adm2_ews$running_z_score_robust
         , typ="l" , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  lines( x = filename_1_adm2_ews$date, y = replicate(length(filename_1_adm2_ews$date),0), col="black",lty=2 )
  legend("topleft",legend=c("adm1","adm2","EWS threshold"),col=c("black","red","black"),lty=c(1,1,2))
  
  #' Difference in robust z-score for leading indicator
  plot( x = filename_1_adm1_ews$date, y = ( filename_1_adm1_ews$running_z_score_robust - filename_1_adm2_ews$running_z_score_robust )
        , xlab = "" , ylab = "z-score difference (adm1 - adm2)"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )
        , ylim = c( max( min( na.omit( filename_1_adm1_ews$running_z_score_robust ) , na.omit( filename_1_adm2_ews$running_z_score_robust ) ) , -10 )  
                    , min( max( na.omit( filename_1_adm1_ews$running_z_score_robust ) , na.omit( filename_1_adm2_ews$running_z_score_robust ) ) , 10 )
                    )
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "blue" #data_colour[2]
  )
  abline(h=0)
  
  #' Repeat for simple LGR mean leading indicator
  #' 
  #' Four plots: hopsitalisations, leading indicator value, z-score, z-score difference
  par( mfrow = c( 4 , 1 ) )
  #' plot hospitalisations
  plot( x = dat_df$date, y = dat_df$cases
        , xlab = "" , ylab = "UK Covid-19 hospitalisations"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )  , ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  #' leading indicator value
  plot( x = filename_2_adm1_ews$date , y = filename_2_adm1_ews$cases
        , xlab = ""  , ylab = "leading indicator value"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) 
        , ylim = c(-0.20, +0.4)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_2_adm2_ews$date , y = filename_2_adm2_ews$cases
         , typ="l"  , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  legend("topleft",legend=c("adm1","adm2"),col=c("black","red"),lty=c(1,1))
  
  #' robust z-score for leading indicator
  plot( x = filename_2_adm1_ews$date  , y = filename_2_adm1_ews$running_z_score_robust
        , xlab = ""  , ylab = "'Robust' z-score of leading indicator"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) 
        , ylim = c( max( min( na.omit( filename_2_adm1_ews$running_z_score_robust ) , na.omit( filename_2_adm2_ews$running_z_score_robust ) ) , -10 )  
                    , min( max( na.omit( filename_2_adm1_ews$running_z_score_robust ) , na.omit( filename_2_adm2_ews$running_z_score_robust ) ) , 10 )
        )
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_2_adm2_ews$date , y = filename_2_adm2_ews$running_z_score_robust
         , typ="l" , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  lines( x = filename_1_adm2_ews$date, y = replicate(length(filename_1_adm2_ews$date),0), col="black",lty=2 )
  legend("topleft",legend=c("adm1","adm2","EWS threshold"),col=c("black","red","black"),lty=c(1,1,2))
  
  #' Difference in robust z-score for leading indicator
  plot( x = filename_2_adm1_ews$date, y = ( filename_2_adm1_ews$running_z_score_robust - filename_2_adm2_ews$running_z_score_robust )
        , xlab = "" , ylab = "z-score difference (adm1 - adm2)"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )
        , ylim = c( max( min( na.omit( filename_2_adm1_ews$running_z_score_robust ) , na.omit( filename_2_adm2_ews$running_z_score_robust ) ) , -10 )  
                    , min( max( na.omit( filename_2_adm1_ews$running_z_score_robust ) , na.omit( filename_2_adm2_ews$running_z_score_robust ) ) , 10 )
        )
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "blue" #data_colour[2]
  )
  abline(h=0)

#### 3 Leading indicator time series  ############################################
  #' 2a: Read in leading indicator time series and calculate 'robust' z-score
  #' Code below adapted from C:\Users\kdrake\OneDrive - Imperial College London\Documents\TFPS\tfps runs\2023_02\quantile_historical\analysis\HPC_scripts\EWS_calc_threshold_lgr_th_085_pval_filter_001_HPC_array.R
  filename_1_adm2_Tg6_5 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm2_Tg6_5/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_001/dataframes_statistics/tfps_pango_lineage_max_lgr.csv"
  filename_1_adm1_Tg6_5 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_001/dataframes_statistics/tfps_pango_lineage_max_lgr.csv"
  filename_2_adm2_Tg6_5 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm2_Tg6_5/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_no/dataframes_statistics/tfps_lgr_simple_mean.csv"
  filename_2_adm1_Tg6_5 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_no/dataframes_statistics/tfps_lgr_simple_mean.csv"
  filename_1_adm2_Tg5 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm2_Tg5/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_001/dataframes_statistics/tfps_pango_lineage_max_lgr.csv"
  filename_1_adm1_Tg5 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm1_Tg5/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_001/dataframes_statistics/tfps_pango_lineage_max_lgr.csv"
  filename_2_adm2_Tg5 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm2_Tg5/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_no/dataframes_statistics/tfps_lgr_simple_mean.csv"
  filename_2_adm1_Tg5 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm1_Tg5/analysis/large_cluster_adjust_lgr_threshold_085/p_val_filter_no/dataframes_statistics/tfps_lgr_simple_mean.csv"
  
  for ( filename in c( filename_1_adm2_Tg6_5 
                      ,filename_1_adm1_Tg6_5 
                      ,filename_2_adm2_Tg6_5 
                      ,filename_2_adm1_Tg6_5 
                      ,filename_1_adm2_Tg5
                      ,filename_1_adm1_Tg5
                      ,filename_2_adm2_Tg5
                      ,filename_2_adm1_Tg5 ) ){
    
    ews_base = data.table::fread( filename ) 
    
    ews_base[ , 1 ] = NULL ; ews_base = data.frame( ews_base )
    
    #' Select relevant column for parameter set required
    if( filename %in% c( filename_1_adm2_Tg6_5 , filename_1_adm1_Tg6_5 , filename_1_adm2_Tg5 , filename_1_adm1_Tg5 ) ){ 
      ews = ews_base[ , c("date","mina07_maxa56_md020")] #' Only need data for one set of parameters for this part
    } else if( filename %in% c( filename_2_adm2_Tg6_5 , filename_2_adm1_Tg6_5 , filename_2_adm2_Tg5 , filename_2_adm1_Tg5 ) ){
      ews = ews_base[ , c("date","mina14_maxa56_md020")] 
    }
    
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
    colnames( ews )[ 5 ] <- c( "running_z_score" )
    colnames( ews )[ 6 ] <- c( "running_z_score_robust" )
    if( filename == filename_1_adm2_Tg6_5 ) { filename_1_adm2_Tg6_5_ews = ews } 
    if( filename == filename_1_adm1_Tg6_5 ) { filename_1_adm1_Tg6_5_ews = ews } 
    if( filename == filename_2_adm2_Tg6_5 ) { filename_2_adm2_Tg6_5_ews = ews } 
    if( filename == filename_2_adm1_Tg6_5 ) { filename_2_adm1_Tg6_5_ews = ews }
    if( filename == filename_1_adm2_Tg5   ) { filename_1_adm2_Tg5_ews = ews } 
    if( filename == filename_1_adm1_Tg5   ) { filename_1_adm1_Tg5_ews = ews } 
    if( filename == filename_2_adm2_Tg5   ) { filename_2_adm2_Tg5_ews = ews } 
    if( filename == filename_2_adm1_Tg5   ) { filename_2_adm1_Tg5_ews = ews } 
  }
  
  
  #' Four plots: hopsitalisations, leading indicator value, z-score, z-score difference
  par( mfrow = c( 4 , 1 ) )
  #' plot hospitalisations
  plot( x = dat_df$date, y = dat_df$cases
        , xlab = "" , ylab = "UK Covid-19 hospitalisations"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )  , ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  #' leading indicator value
  plot( x = filename_1_adm1_Tg6_5_ews$date , y = filename_1_adm1_Tg6_5_ews$cases
        , xlab = ""  , ylab = "leading indicator value"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_1_adm1_Tg5_ews$date , y = filename_1_adm1_Tg5_ews$cases
         , xlab = ""  , ylab = "leading indicator value"
         , cex.lab = 1.25  , cex.axis = 1.25
         , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
         #, xaxt = "n"
         , xaxs="i" , yaxs="i" #' Makes plot box tight to values
         , typ="l"  , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  legend("topleft",legend=c("Tg=6.5 days (adm1)","Tg=5 days (adm1)"),col=c("black","red"),lty=c(1,1))
  
  #' robust z-score for leading indicator
  plot( x = filename_1_adm1_Tg6_5_ews$date  , y = filename_1_adm1_Tg6_5_ews$running_z_score_robust
        , xlab = ""  , ylab = "'Robust' z-score of leading indicator"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_1_adm1_Tg5_ews$date , y = filename_1_adm1_Tg5_ews$running_z_score_robust
         , typ="l" , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  lines( x = filename_1_adm1_Tg6_5_ews$date, y = replicate(length(filename_1_adm1_Tg6_5_ews$date),0), col="black",lty=2 )
  legend("topleft",legend=c("Tg=6.5 days (adm1)","Tg=5 days (adm1)","EWS threshold"),col=c("black","red","black"),lty=c(1,1,2))
  
  #' Difference in robust z-score for leading indicator
  plot( x = filename_1_adm1_Tg6_5_ews$date, y = ( filename_1_adm1_Tg6_5_ews$running_z_score_robust - filename_1_adm1_Tg5_ews$running_z_score_robust )
        , xlab = "" , ylab = "z-score difference (adm1 - adm2)"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )
        , ylim = c( max( min( na.omit( filename_1_adm1_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_1_adm1_Tg5_ews$running_z_score_robust ) ) , -10 )  
                    , min( max( na.omit( filename_1_adm1_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_1_adm1_Tg5_ews$running_z_score_robust ) ) , 10 )
        )
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "blue" #data_colour[2]
  )
  abline(h=0)
  
  #' Repeat for simple LGR mean leading indicator
  #' 
  #' Four plots: hopsitalisations, leading indicator value, z-score, z-score difference
  par( mfrow = c( 4 , 1 ) )
  #' plot hospitalisations
  plot( x = dat_df$date, y = dat_df$cases
        , xlab = "" , ylab = "UK Covid-19 hospitalisations"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )  , ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  #' leading indicator value
  plot( x = filename_2_adm1_Tg6_5_ews$date , y = filename_2_adm1_Tg6_5_ews$cases
        , xlab = ""  , ylab = "leading indicator value"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_2_adm1_Tg5_ews$date , y = filename_2_adm1_Tg5_ews$cases
         , xlab = ""  , ylab = "leading indicator value"
         , cex.lab = 1.25  , cex.axis = 1.25
         , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
         #, xaxt = "n"
         , xaxs="i" , yaxs="i" #' Makes plot box tight to values
         , typ="l"  , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  legend("bottomleft",legend=c("Tg=6.5 days (adm1)","Tg=5 days (adm1)"),col=c("black","red"),lty=c(1,1))
  
  #' robust z-score for leading indicator
  plot( x = filename_2_adm1_Tg6_5_ews$date  , y = filename_2_adm1_Tg6_5_ews$running_z_score_robust
        , xlab = ""  , ylab = "'Robust' z-score of leading indicator"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_2_adm1_Tg5_ews$date , y = filename_2_adm1_Tg5_ews$running_z_score_robust
         , typ="l" , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  lines( x = filename_2_adm1_Tg6_5_ews$date, y = replicate(length(filename_2_adm1_Tg6_5_ews$date),0), col="black",lty=2 )
  legend("bottomleft",legend=c("Tg=6.5 days (adm1)","Tg=5 days (adm1)","EWS threshold"),col=c("black","red","black"),lty=c(1,1,2))
  
  #' Difference in robust z-score for leading indicator
  plot( x = filename_2_adm1_Tg6_5_ews$date, y = ( filename_2_adm1_Tg6_5_ews$running_z_score_robust - filename_2_adm1_Tg5_ews$running_z_score_robust )
        , xlab = "" , ylab = "z-score difference (adm1 - adm2)"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )
        , ylim = c( max( min( na.omit( filename_2_adm1_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_2_adm1_Tg5_ews$running_z_score_robust ) ) , -10 )  
                    , min( max( na.omit( filename_2_adm1_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_2_adm1_Tg5_ews$running_z_score_robust ) ) , 10 )
        )
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "blue" #data_colour[2]
  )
  abline(h=0)
  
  #' Repeat for adm2 
  #' Pango lineage max LGR leading indicator
  #' 
  #' Four plots: hopsitalisations, leading indicator value, z-score, z-score difference
  par( mfrow = c( 4 , 1 ) )
  #' plot hospitalisations
  plot( x = dat_df$date, y = dat_df$cases
        , xlab = "" , ylab = "UK Covid-19 hospitalisations"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )  , ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  #' leading indicator value
  plot( x = filename_1_adm2_Tg6_5_ews$date , y = filename_1_adm2_Tg6_5_ews$cases
        , xlab = ""  , ylab = "leading indicator value"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_1_adm2_Tg5_ews$date , y = filename_1_adm2_Tg5_ews$cases
         , xlab = ""  , ylab = "leading indicator value"
         , cex.lab = 1.25  , cex.axis = 1.25
         , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
         #, xaxt = "n"
         , xaxs="i" , yaxs="i" #' Makes plot box tight to values
         , typ="l"  , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  legend("topleft",legend=c("Tg=6.5 days (adm2)","Tg=5 days (adm2)"),col=c("black","red"),lty=c(1,1))
  
  #' robust z-score for leading indicator
  plot( x = filename_1_adm2_Tg6_5_ews$date  , y = filename_1_adm2_Tg6_5_ews$running_z_score_robust
        , xlab = ""  , ylab = "'Robust' z-score of leading indicator"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) 
        , ylim = c( max( min( na.omit( filename_1_adm2_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_1_adm2_Tg5_ews$running_z_score_robust ) ) , -10 )  
                    , min( max( na.omit( filename_1_adm2_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_1_adm2_Tg5_ews$running_z_score_robust ) ) , 10 ) )
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_1_adm2_Tg5_ews$date , y = filename_1_adm2_Tg5_ews$running_z_score_robust
         , typ="l" , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  lines( x = filename_1_adm2_Tg6_5_ews$date, y = replicate(length(filename_1_adm2_Tg6_5_ews$date),0), col="black",lty=2 )
  legend("bottomleft",legend=c("Tg=6.5 days (adm2)","Tg=5 days (adm2)","EWS threshold"),col=c("black","red","black"),lty=c(1,1,2))
  
  #' Difference in robust z-score for leading indicator
  plot( x = filename_1_adm2_Tg6_5_ews$date, y = ( filename_1_adm2_Tg6_5_ews$running_z_score_robust - filename_1_adm2_Tg5_ews$running_z_score_robust )
        , xlab = "" , ylab = "z-score difference (Tg6.5 - Tg5 (adm2))"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )
        , ylim = c( max( min( na.omit( filename_1_adm2_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_1_adm2_Tg5_ews$running_z_score_robust ) ) , -10 )  
                    , min( max( na.omit( filename_1_adm2_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_1_adm2_Tg5_ews$running_z_score_robust ) ) , 10 )
        )
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "blue" #data_colour[2]
  )
  abline(h=0)
  
  #' Repeat for adm2 
  #' Repeat for simple LGR mean leading indicator
  #' 
  #' Four plots: hopsitalisations, leading indicator value, z-score, z-score difference
  par( mfrow = c( 4 , 1 ) )
  #' plot hospitalisations
  plot( x = dat_df$date, y = dat_df$cases
        , xlab = "" , ylab = "UK Covid-19 hospitalisations"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )  , ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  #' leading indicator value
  plot( x = filename_2_adm2_Tg6_5_ews$date , y = filename_2_adm2_Tg6_5_ews$cases
        , xlab = ""  , ylab = "leading indicator value"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_2_adm2_Tg5_ews$date , y = filename_2_adm2_Tg5_ews$cases
         , xlab = ""  , ylab = "leading indicator value"
         , cex.lab = 1.25  , cex.axis = 1.25
         , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) #, ylim = c(0, 5000)
         #, xaxt = "n"
         , xaxs="i" , yaxs="i" #' Makes plot box tight to values
         , typ="l"  , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  legend("bottomleft",legend=c("Tg=6.5 days (adm2)","Tg=5 days (adm2)"),col=c("black","red"),lty=c(1,1))
  
  #' robust z-score for leading indicator
  plot( x = filename_2_adm2_Tg6_5_ews$date  , y = filename_2_adm2_Tg6_5_ews$running_z_score_robust
        , xlab = ""  , ylab = "'Robust' z-score of leading indicator"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) ) 
        , ylim = c( max( min( na.omit( filename_2_adm2_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_2_adm2_Tg5_ews$running_z_score_robust ) ) , -10 )  
                    , min( max( na.omit( filename_2_adm2_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_2_adm2_Tg5_ews$running_z_score_robust ) ) , 10 ) )
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l" , lwd = 1.5
        , col = "black" #data_colour[2]
  )
  lines( x = filename_2_adm2_Tg5_ews$date , y = filename_2_adm2_Tg5_ews$running_z_score_robust
         , typ="l" , lwd = 1.5
         , col = "red" #data_colour[2]
  )
  lines( x = filename_2_adm2_Tg6_5_ews$date, y = replicate(length(filename_2_adm2_Tg6_5_ews$date),0), col="black",lty=2 )
  legend("bottomleft",legend=c("Tg=6.5 days (adm2)","Tg=5 days (adm2)","EWS threshold"),col=c("black","red","black"),lty=c(1,1,2))
  
  #' Difference in robust z-score for leading indicator
  plot( x = filename_2_adm2_Tg6_5_ews$date, y = ( filename_2_adm2_Tg6_5_ews$running_z_score_robust - filename_2_adm2_Tg5_ews$running_z_score_robust )
        , xlab = "" , ylab = "z-score difference (Tg6.5 - Tg5 (adm2))"
        , cex.lab = 1.25 , cex.axis = 1.25
        , xlim = c( dat_df$date[ 1 ] , max( dat_df$date ) )
        , ylim = c( max( min( na.omit( filename_2_adm2_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_2_adm2_Tg5_ews$running_z_score_robust ) ) , -10 )  
                    , min( max( na.omit( filename_2_adm2_Tg6_5_ews$running_z_score_robust ) , na.omit( filename_2_adm2_Tg5_ews$running_z_score_robust ) ) , 10 )
        )
        #, xaxt = "n"
        , xaxs="i" , yaxs="i" #' Makes plot box tight to values
        , typ="l"  , lwd = 1.5
        , col = "blue" #data_colour[2]
  )
  abline(h=0)
  
#### 4: Extract data on lead times and false positives ########################

#' Read in dataframes

#' Analysis using adm2 region and generation time (Tg) = 6.5 days
#' This file looks like it is filtered for parameters that produce an EWS for all waves
wave_results_analysis_reshaped_fp_df_adm2_Tg6_5 = readRDS("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm2_Tg6_5/analysis/EWS_df_combined/wave_results_analysis_reshaped_fp_df.rds")
#' and this one doesn't (also different formats)
wave_results_analysis_reshaped_fp_adm2_Tg6_5 = readRDS("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm2_Tg6_5/analysis/EWS_df_combined/wave_results_analysis_reshaped_fp.rds")
  
#'format dates
wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w2_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w2_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w3_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w3_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w4_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w4_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w5_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w5_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w6_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w6_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w7_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w7_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w8_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$w8_earliest_tp_EWS_date , origin = "1970-01-01" )

wave_results_analysis_reshaped_fp_adm2_Tg6_5$w2_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg6_5$w2_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg6_5$w3_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg6_5$w3_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg6_5$w4_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg6_5$w4_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg6_5$w5_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg6_5$w5_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg6_5$w6_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg6_5$w6_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg6_5$w7_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg6_5$w7_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg6_5$w8_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg6_5$w8_earliest_tp_EWS_date , origin = "1970-01-01" )

#' subset for leading indicators / parameters of interest
wrar_fp_df_adm2_Tg6_5 = as.data.frame( t( subset( wave_results_analysis_reshaped_fp_df_adm2_Tg6_5 
                                , wave_results_analysis_reshaped_fp_df_adm2_Tg6_5$leading_indicator_type %in% c("tfps_pango_lineage_max_lgr" , "tfps_lgr_simple_mean")
                                )))
View(wrar_fp_df_adm2_Tg6_5)
wrar_fp_df_adm2_Tg6_5 = wrar_fp_df_adm2_Tg6_5[,c(3,4)]
colnames(wrar_fp_df_adm2_Tg6_5)

wrar_fp_adm2_Tg6_5 = as.data.frame( t( subset( wave_results_analysis_reshaped_fp_adm2_Tg6_5 
                                                  , wave_results_analysis_reshaped_fp_adm2_Tg6_5$variables_code %in% c("tfps_pango_lineage_max_lgr_07_56_20_0.01_0.85_0" , "tfps_lgr_simple_mean_14_56_20_10000_0.85_0")
)))
View(wrar_fp_adm2_Tg6_5)
#wrar_fp_adm2_Tg6_5 = wrar_fp_adm2_Tg6_5[,c(3,4)]
#colnames(wrar_fp_df_adm2_Tg6_5)


#' Analysis using adm1 region and generation time (Tg) = 6.5 days
#' This file looks like it is filtered for parameters that produce an EWS for all waves
wave_results_analysis_reshaped_fp_df_adm1_Tg6_5 = readRDS("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined_2023_03_20/wave_results_analysis_reshaped_fp_df.rds")
#' and this one doesn't (also different formats)
wave_results_analysis_reshaped_fp_adm1_Tg6_5 = readRDS("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined_2023_03_20/wave_results_analysis_reshaped_fp.rds")

#'format dates
wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w2_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w2_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w3_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w3_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w4_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w4_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w5_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w5_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w6_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w6_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w7_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w7_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w8_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$w8_earliest_tp_EWS_date , origin = "1970-01-01" )
#' subset for leading indicators / parameters of interest
wrar_fp_df_adm1_Tg6_5 = as.data.frame( t( subset( wave_results_analysis_reshaped_fp_df_adm1_Tg6_5 
                                                  , wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$leading_indicator_type %in% c("tfps_pango_lineage_max_lgr" , "tfps_lgr_simple_mean") &
                                                    wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$TFPS_cluster_min_age %in% c("07" , "14") &
                                                    wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$EWS_threshold == 0.00 &
                                                    wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$TFPS_cluster_max_age %in% c("56") &
                                                    wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$TFPS_cluster_min_descendants %in% c("20") &
                                                    wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$parent_sub_lgr_threshold == 0.85 &
                                                    wave_results_analysis_reshaped_fp_df_adm1_Tg6_5$LGR_p_value_threshold %in% c(10000,0.01)
)))
View( wrar_fp_df_adm1_Tg6_5 )
wrar_fp_df_adm1_Tg6_5 = wrar_fp_df_adm1_Tg6_5[,c(4,5)]
colnames( wrar_fp_df_adm1_Tg6_5 )

#' Analysis using adm2 region and generation time (Tg) = 5 days
#' This file looks like it is filtered for parameters that produce an EWS for all waves
wave_results_analysis_reshaped_fp_df_adm2_Tg5 = readRDS("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm2_Tg5/analysis/EWS_df_combined/wave_results_analysis_reshaped_fp_df.rds")
#' and this one doesn't (also different formats)
wave_results_analysis_reshaped_fp_adm2_Tg5 = readRDS("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm2_Tg5/analysis/EWS_df_combined/wave_results_analysis_reshaped_fp.rds")

#'format dates
wave_results_analysis_reshaped_fp_df_adm2_Tg5$w2_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg5$w2_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg5$w3_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg5$w3_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg5$w4_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg5$w4_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg5$w5_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg5$w5_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg5$w6_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg5$w6_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg5$w7_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg5$w7_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm2_Tg5$w8_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm2_Tg5$w8_earliest_tp_EWS_date , origin = "1970-01-01" )

wave_results_analysis_reshaped_fp_adm2_Tg5$w2_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg5$w2_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg5$w3_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg5$w3_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg5$w4_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg5$w4_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg5$w5_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg5$w5_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg5$w6_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg5$w6_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg5$w7_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg5$w7_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_adm2_Tg5$w8_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_adm2_Tg5$w8_earliest_tp_EWS_date , origin = "1970-01-01" )

#' subset for leading indicators / parameters of interest
wrar_fp_df_adm2_Tg5 = as.data.frame( t( subset( wave_results_analysis_reshaped_fp_df_adm2_Tg5 
                                                  , wave_results_analysis_reshaped_fp_df_adm2_Tg5$leading_indicator_type %in% c("tfps_pango_lineage_max_lgr" , "tfps_lgr_simple_mean") #&
                                                    #wave_results_analysis_reshaped_fp_df_adm2_Tg5$TFPS_cluster_min_age %in% c("07" , "14") &
                                                    #wave_results_analysis_reshaped_fp_df_adm2_Tg5$EWS_threshold == 0.00 &
                                                    #wave_results_analysis_reshaped_fp_df_adm2_Tg5$TFPS_cluster_max_age %in% c("56") &
                                                    #wave_results_analysis_reshaped_fp_df_adm2_Tg5$TFPS_cluster_min_descendants %in% c("20") &
                                                    #wave_results_analysis_reshaped_fp_df_adm2_Tg5$parent_sub_lgr_threshold == 0.85 &
                                                    #wave_results_analysis_reshaped_fp_df_adm2_Tg5$LGR_p_value_threshold %in% c(1000,0.01)
)))
View( wrar_fp_df_adm2_Tg5 )
#wrar_fp_df_adm2_Tg5 = wrar_fp_df_adm2_Tg5[,c(2,3)]
colnames( wrar_fp_df_adm2_Tg5 )

wrar_fp_adm2_Tg5 = as.data.frame( t( subset( wave_results_analysis_reshaped_fp_adm2_Tg5 
                                               , wave_results_analysis_reshaped_fp_adm2_Tg5$variables_code %in% c("tfps_pango_lineage_max_lgr_07_56_20_0.01_0.85_0" , "tfps_lgr_simple_mean_14_56_20_10000_0.85_0")
)))
View(wrar_fp_adm2_Tg5)


#' Analysis using adm1 region and generation time (Tg) = 5 days
wave_results_analysis_reshaped_fp_df_adm1_Tg5 = readRDS("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/adm1_Tg5/analysis/EWS_df_combined/wave_results_analysis_reshaped_fp_df.rds")

#'format dates
wave_results_analysis_reshaped_fp_df_adm1_Tg5$w2_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg5$w2_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg5$w3_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg5$w3_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg5$w4_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg5$w4_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg5$w5_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg5$w5_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg5$w6_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg5$w6_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg5$w7_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg5$w7_earliest_tp_EWS_date , origin = "1970-01-01" )
wave_results_analysis_reshaped_fp_df_adm1_Tg5$w8_earliest_tp_EWS_date = as.Date( wave_results_analysis_reshaped_fp_df_adm1_Tg5$w8_earliest_tp_EWS_date , origin = "1970-01-01" )
#' subset for leading indicators / parameters of interest
wrar_fp_df_adm1_Tg5 = as.data.frame( t( subset( wave_results_analysis_reshaped_fp_df_adm1_Tg5 
                                                , wave_results_analysis_reshaped_fp_df_adm1_Tg5$leading_indicator_type %in% c("tfps_pango_lineage_max_lgr" , "tfps_lgr_simple_mean") #&
                                                #wave_results_analysis_reshaped_fp_df_adm2_Tg5$TFPS_cluster_min_age %in% c("07" , "14") &
                                                #wave_results_analysis_reshaped_fp_df_adm2_Tg5$EWS_threshold == 0.00 &
                                                #wave_results_analysis_reshaped_fp_df_adm2_Tg5$TFPS_cluster_max_age %in% c("56") &
                                                #wave_results_analysis_reshaped_fp_df_adm2_Tg5$TFPS_cluster_min_descendants %in% c("20") &
                                                #wave_results_analysis_reshaped_fp_df_adm2_Tg5$parent_sub_lgr_threshold == 0.85 &
                                                #wave_results_analysis_reshaped_fp_df_adm2_Tg5$LGR_p_value_threshold %in% c(1000,0.01)
)))
View( wrar_fp_df_adm1_Tg5 )
wrar_fp_df_adm1_Tg5 = wrar_fp_df_adm1_Tg5[,c(4,5)]
colnames( wrar_fp_df_adm1_Tg5 )


#' Comparison of all leading indicator / parameter sets run using adm2 against adm1 run
res_adm2 = wave_results_analysis_reshaped_fp_adm2_Tg6_5
res_adm1 = subset( wave_results_analysis_reshaped_fp_adm1_Tg6_5 , wave_results_analysis_reshaped_fp_adm1_Tg6_5$variables_code %in% res_adm2$variables_code ) #' filter for same parameter sets run using adm2
res_adm2 = res_adm2[2:nrow(res_adm2),] #' Remove 1st row as it is NA
res_adm1 = res_adm1[2:nrow(res_adm1),] #' Remove 1st row as it is NA
#' Plot histograms
par(mfrow=c(1,1))
#' wave 2 B.1.177
hist(res_adm1$w2_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w2_lead_lag_days),na.omit(res_adm2$w2_lead_lag_days)) , max(na.omit(res_adm1$w2_lead_lag_days),na.omit(res_adm2$w2_lead_lag_days)) ) , ylim=c(0,50) , col=rgb(0,0,1,1/4) , main = "B.1.177" , xlab = "EWS lead (-ve) or lag(+ve) time (days)" )
hist(res_adm2$w2_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w2_lead_lag_days),na.omit(res_adm2$w2_lead_lag_days)) , max(na.omit(res_adm1$w2_lead_lag_days),na.omit(res_adm2$w2_lead_lag_days)) ) , ylim = c(0,50), col=rgb(1,0,0,1/4) , add=T)
legend("topright",legend=c("adm1","adm2"),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), pch=c(15,15))
#' wave 3 alpha
hist(res_adm1$w3_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w3_lead_lag_days),na.omit(res_adm2$w3_lead_lag_days)) , max(na.omit(res_adm1$w3_lead_lag_days),na.omit(res_adm2$w3_lead_lag_days)) ) , ylim=c(0,30) , col=rgb(0,0,1,1/4) , main = "Alpha" , xlab = "EWS lead (-ve) or lag(+ve) time (days)" )
hist(res_adm2$w3_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w3_lead_lag_days),na.omit(res_adm2$w3_lead_lag_days)) , max(na.omit(res_adm1$w3_lead_lag_days),na.omit(res_adm2$w3_lead_lag_days)) ) , ylim = c(0,30), col=rgb(1,0,0,1/4) , add=T)
legend("topright",legend=c("adm1","adm2"),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), pch=c(15,15))
#' wave 4 Delta (1)
hist(res_adm1$w4_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w4_lead_lag_days),na.omit(res_adm2$w4_lead_lag_days)) , max(na.omit(res_adm1$w4_lead_lag_days),na.omit(res_adm2$w4_lead_lag_days)) ) , ylim=c(0,100) , col=rgb(0,0,1,1/4) , main = "Delta (1)" , xlab = "EWS lead (-ve) or lag(+ve) time (days)" )
hist(res_adm2$w4_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w4_lead_lag_days),na.omit(res_adm2$w4_lead_lag_days)) , max(na.omit(res_adm1$w4_lead_lag_days),na.omit(res_adm2$w4_lead_lag_days)) ) , ylim = c(0,100), col=rgb(1,0,0,1/4) , add=T)
legend("topright",legend=c("adm1","adm2"),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), pch=c(15,15))
#' wave 5 Delta (2)
hist(res_adm1$w5_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w5_lead_lag_days),na.omit(res_adm2$w5_lead_lag_days)) , max(na.omit(res_adm1$w5_lead_lag_days),na.omit(res_adm2$w5_lead_lag_days)) ) , ylim=c(0,100) , col=rgb(0,0,1,1/4) , main = "Delta (2)" , xlab = "EWS lead (-ve) or lag(+ve) time (days)" )
hist(res_adm2$w5_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w5_lead_lag_days),na.omit(res_adm2$w5_lead_lag_days)) , max(na.omit(res_adm1$w5_lead_lag_days),na.omit(res_adm2$w5_lead_lag_days)) ) , ylim = c(0,100), col=rgb(1,0,0,1/4) , add=T)
legend("topright",legend=c("adm1","adm2"),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), pch=c(15,15))
#' wave 6 Delta (3)
hist(res_adm1$w6_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w6_lead_lag_days),na.omit(res_adm2$w6_lead_lag_days)) , max(na.omit(res_adm1$w6_lead_lag_days),na.omit(res_adm2$w6_lead_lag_days)) ) , ylim=c(0,40) , col=rgb(0,0,1,1/4) , main = "Delta (3)" , xlab = "EWS lead (-ve) or lag(+ve) time (days)" )
hist(res_adm2$w6_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w6_lead_lag_days),na.omit(res_adm2$w6_lead_lag_days)) , max(na.omit(res_adm1$w6_lead_lag_days),na.omit(res_adm2$w6_lead_lag_days)) ) , ylim = c(0,40), col=rgb(1,0,0,1/4) , add=T)
legend("topright",legend=c("adm1","adm2"),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), pch=c(15,15))
#' wave 7 Omicron BA.1
hist(res_adm1$w7_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w7_lead_lag_days),na.omit(res_adm2$w7_lead_lag_days)) , max(na.omit(res_adm1$w7_lead_lag_days),na.omit(res_adm2$w7_lead_lag_days)) ) , ylim=c(0,50) , col=rgb(0,0,1,1/4) , main = "Omicron BA.1" , xlab = "EWS lead (-ve) or lag(+ve) time (days)" )
hist(res_adm2$w7_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w7_lead_lag_days),na.omit(res_adm2$w7_lead_lag_days)) , max(na.omit(res_adm1$w7_lead_lag_days),na.omit(res_adm2$w7_lead_lag_days)) ) , ylim = c(0,50), col=rgb(1,0,0,1/4) , add=T)
legend("topright",legend=c("adm1","adm2"),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), pch=c(15,15))
#' wave 8 Omicron BA.2
hist(res_adm1$w8_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w8_lead_lag_days),na.omit(res_adm2$w8_lead_lag_days)) , max(na.omit(res_adm1$w8_lead_lag_days),na.omit(res_adm2$w8_lead_lag_days)) ) , ylim=c(0,50) , col=rgb(0,0,1,1/4) , main = "Omicron BA.2" , xlab = "EWS lead (-ve) or lag(+ve) time (days)" )
hist(res_adm2$w8_lead_lag_days,breaks = 100, xlim = c( min(na.omit(res_adm1$w8_lead_lag_days),na.omit(res_adm2$w8_lead_lag_days)) , max(na.omit(res_adm1$w8_lead_lag_days),na.omit(res_adm2$w8_lead_lag_days)) ) , ylim = c(0,50), col=rgb(1,0,0,1/4) , add=T)
legend("topright",legend=c("adm1","adm2"),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), pch=c(15,15))

#' mean across waves
lead_times_mean = data.frame("test")
lead_times_mean$w2_adm1 = mean(res_adm1$w2_lead_lag_days) # +10
lead_times_mean$w2_adm2 = mean(res_adm2$w2_lead_lag_days) # +7
lead_times_mean$X.test. = NULL

lead_times_mean$w3_adm1 = mean( na.omit(res_adm1$w3_lead_lag_days) ) # +21.1
table(is.na( res_adm1$w3_lead_lag_days ))[2] # 18 NAs
lead_times_mean$w3_adm2 = mean( na.omit(res_adm2$w3_lead_lag_days) ) # +14.9
table(is.na( res_adm2$w3_lead_lag_days ))[2] # 44 NAs

lead_times_mean$w4_adm1 = mean( na.omit(res_adm1$w4_lead_lag_days) ) # -18.6
table(is.na( res_adm1$w4_lead_lag_days ))[2] # 0 NAs
lead_times_mean$w4_adm2 = mean( na.omit(res_adm2$w4_lead_lag_days) ) # -21.2
table(is.na( res_adm2$w4_lead_lag_days ))[2] # 0 NAs

lead_times_mean$w5_adm1 = mean( na.omit(res_adm1$w5_lead_lag_days) ) # +0.2
table(is.na( res_adm1$w5_lead_lag_days ))[2] # 0 NAs
lead_times_mean$w5_adm2 = mean( na.omit(res_adm2$w5_lead_lag_days) ) # +1.1
table(is.na( res_adm2$w5_lead_lag_days ))[2] # 7 NAs

lead_times_mean$w6_adm1 = mean( na.omit(res_adm1$w6_lead_lag_days) ) # -1.6
table(is.na( res_adm1$w6_lead_lag_days ))[2] # 4 NAs
lead_times_mean$w6_adm2 = mean( na.omit(res_adm2$w6_lead_lag_days) ) # -2.6
table(is.na( res_adm2$w6_lead_lag_days ))[2] # 4 NAs

lead_times_mean$w7_adm1 = mean( na.omit(res_adm1$w7_lead_lag_days) ) # +16.8
table(is.na( res_adm1$w7_lead_lag_days ))[2] # 0 NAs
lead_times_mean$w7_adm2 = mean( na.omit(res_adm2$w7_lead_lag_days) ) # +16.8
table(is.na( res_adm2$w7_lead_lag_days ))[2] # 0 NAs

lead_times_mean$w8_adm1 = mean( na.omit(res_adm1$w8_lead_lag_days) ) # -12.3
table(is.na( res_adm1$w8_lead_lag_days ))[2] # 0 NAs
lead_times_mean$w8_adm2 = mean( na.omit(res_adm2$w8_lead_lag_days) ) # -12.3
table(is.na( res_adm2$w8_lead_lag_days ))[2] # 0 NAs

#' Compute mean improvement / worsening per wave and in total by switching from adm1 to adm2
res_merge = merge( x = res_adm2 , y = res_adm1
                   , by.x = "variables_code" , by.y = "variables_code" 
                   , all.x = TRUE
                   , all.y = TRUE
                   )
res_merge$w2_lead_lag_change = res_merge$w2_lead_lag_days.y - res_merge$w2_lead_lag_days.x #' lead or lag using adm1 - lead or lag using adm2
res_merge$w3_lead_lag_change = res_merge$w3_lead_lag_days.y - res_merge$w3_lead_lag_days.x #' lead or lag using adm1 - lead or lag using adm2
res_merge$w4_lead_lag_change = res_merge$w4_lead_lag_days.y - res_merge$w4_lead_lag_days.x #' lead or lag using adm1 - lead or lag using adm2
res_merge$w5_lead_lag_change = res_merge$w5_lead_lag_days.y - res_merge$w5_lead_lag_days.x #' lead or lag using adm1 - lead or lag using adm2
res_merge$w6_lead_lag_change = res_merge$w6_lead_lag_days.y - res_merge$w6_lead_lag_days.x #' lead or lag using adm1 - lead or lag using adm2
res_merge$w7_lead_lag_change = res_merge$w7_lead_lag_days.y - res_merge$w7_lead_lag_days.x #' lead or lag using adm1 - lead or lag using adm2
res_merge$w8_lead_lag_change = res_merge$w8_lead_lag_days.y - res_merge$w8_lead_lag_days.x #' lead or lag using adm1 - lead or lag using adm2

#### 5: Save environment - used in plots for EWS paper ########
save.image( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/comparison_adm_Tg.RData" )
