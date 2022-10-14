#' Author: Kieran Drake
#' Date: May 2022
#' Purpose: Calculation of early warning signal (EWS)
#'  Based on methods described in 
#'  O’Brien, D. A., & Clements, C. F. (2021).
#'  Early warning signal reliability varies with COVID-19 waves
#'  https://doi.org/10.6084/M9.FIGSHARE.17081673.V1

################## Load libraries ##############################################
# install packages for fread (which is quicker than read_csv)
if(!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") 
if(!requireNamespace("Rtools", quietly = TRUE))
  install.packages("Rtools")
if(!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
if(!requireNamespace("R.utils", quietly = TRUE))
  install.packages('R.utils')
#if(!requireNamespace("TTR", quietly = TRUE))
#install.packages("TTR")
#if(!requireNamespace("gratia", quietly = TRUE))
#install.packages("gratia") #https://gavinsimpson.github.io/gratia/

# Packages to calculate skewness
if(!requireNamespace("e1071", quietly = TRUE))
  install.packages('e1071')
if(!requireNamespace("asbio", quietly = TRUE))
  install.packages('asbio')
if(!requireNamespace("prettyR", quietly = TRUE))
  install.packages('prettyR')

library(BiocManager)
library(data.table)
library(R.utils)
library(e1071)
library(asbio)
library(prettyR)

################## Load Functions #########################################
#' Functions to load
#' data_load() @ C:\Users\kdrake\OneDrive - Imperial College London\Documents\GitHub\Early_Warning_Signal
#' GAM_fitting() @ C:\Users\kdrake\OneDrive - Imperial College London\Documents\GitHub\Early_Warning_Signal
#' wave_reset_derivative_method @ C:\Users\kdrake\OneDrive - Imperial College London\Documents\GitHub\Early_Warning_Signal

################## Load data ##################################################
#' Load case/hospitalisation data
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs'
#filename <- "data_2022-May-09 - UK Covid-19 cases by sample date.csv"
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
#dat_type <- "cases"
dat_type <- "hospitalisations"
#' data_load() @ C:\Users\kdrake\GitHub\Early_Warning_Signal
dat_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)

#' Load leading indicators
#' Select leading indicator
#' 1 = Test data - shifted hospitalisations
#' 2 = variance of cluster logistic growth rates
#' 3 = SARS-CoV-2 Ct values (PCR cycle threshold)
#' 4 = SARS-CoV-2 positivity rates
#' 5 = Behavioural - CoMix survey
#' 6 = Behavioural - Google mobility
lead_ind_type = 2 #

#' 1 - Test data - shifted hospitalisations
if ( lead_ind_type == 1 ){
  shift <- function( x , n ){
    c( x[ -( seq( n ) ) ], rep( NA , n ) )
  }
  
  ews_base = data.frame(   "date" = dat_df$date 
                         , "hosp_shift_00" = dat_df$cases
                         , "hosp_shift_10" = shift(dat_df$cases, 10)
                         , "hosp_shift_20" = shift(dat_df$cases, 20) 
                         , "hosp_shift_30" = shift(dat_df$cases, 30) 
                         )
  filename_prefix = "test_hosp_shift_"
}
#' 2 - variance of cluster logistic growth rates
if ( lead_ind_type == 2 ){
  folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/analysis/p_val_filter_no/dataframes - statistics'
  filenames <- c(   "tfps_vlgr_samp.csv"    , "tfps_vlgr_simple_samp.csv"    , "tfps_vlgr_gam_samp.csv" 
                  , "tfps_vlgr_pop.csv"     , "tfps_vlgr_simple_pop.csv"     , "tfps_vlgr_gam_pop.csv"
                  , "tfps_vlgr_wtd.csv"
                  , "tfps_lgr_max.csv"      , "tfps_lgr_simple_max.csv"      , "tfps_lgr_gam_max.csv"
                  , "tfps_lgr_mean.csv"     , "tfps_lgr_simple_mean.csv"     , "tfps_lgr_gam_mean.csv"
                  , "tfps_lgr_wtd_mean.csv" , "tfps_lgr_simple_wtd_mean.csv" , "tfps_lgr_gam_wtd_mean.csv"
                  , "tfps_clock_outlier_max.csv" , "tfps_clock_outlier_mean.csv"
                  )
  filenames_prefix <- paste( data.frame( str_split( filenames , pattern = ".csv" ) )[ 1 , ] , "_" , sep = "" )
  
  #filename <- "tfps_vlgr_mdperc.csv" #"tfps_vlgr_mdperc.csv" #"tfps_lead_ind_comp.csv" # "tfps_vlgr.csv" "tfps_v_gam_lgr.csv" "tfps_vlgr_wtd.csv" "tfps_lgr_max.csv"
  #filename_prefix = "vlgr_mdperc_" # "vlgr_md_perc" #"lead_ind_comp_" #"lgr_max_"
  setwd( folder )
  #ews_base = fread( filename ) ;  ews_base[,1] = NULL ; ews_base = data.frame( ews_base )
}
#' 3 - SARS-CoV-2 Ct values (PCR cycle threshold)
if ( lead_ind_type == 3 ){
  folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files'
  #filename <- "Ct_p2_mean_df.csv"
  #filename_prefix = "Ct_p2_mean_"
  #' Or 
  filename <- "Ct_p2_median_df.csv"
  filename_prefix = "Ct_p2_median_"
  setwd( folder )
  ews_base = fread( filename ) ;  ews_base[,1] = NULL ; ews_base = data.frame( ews_base )
  # Select data to use for EWS calculation
  #ews_y <- c("O_Ct", "N_Ct", "S_Ct", "Control_Ct", "Ct_min", "Ct_mean", "O_Ct_norm", "N_Ct_norm", "S_Ct_norm", "O_vl", "N_vl", "S_vl", "vl_min", "vl_mean", "Ct_min_skew", "Ct_min_stdev", "vl_max_skew", "vl_max_stdev")
  #ews_x <- "Date" 
  #ews = ews[,c("Date","Ct_min_stdev")] #O_Ct, N_Ct, S_Ct, Control_Ct, Ct_min, Ct_mean, O_Ct_norm, N_Ct_norm, S_Ct_norm, O_vl, N_vl, S_vl, vl_min, vl_mean, Ct_min_skew, Ct_min_stdev, vl_max_skew, vl_max_stdev
  #dat_type = "O-gene Ct values"
}
#' 4 - SARS-CoV-2 positivity rates
if ( lead_ind_type == 4 ){
  folder <- "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/Positivity rates/Outputs/"
  filename <- "sc2_positivity_rate_England.rds"
  filename_prefix = "PCR_positivity_rate_"
  setwd( folder )
  ews_base = readRDS( filename ) ; ews_base = data.frame( ews_base[,c(1,4,5)] )
  
}
#' 5 - Behavioural - CoMix survey
if ( lead_ind_type == 5 ){
  folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/UK CoMix data 2022/CoMix data'
  filename <- "2022-03-02_bs_means_2w_open.csv"
  filename_prefix = "CoMix_" 
  setwd( folder )
  ews_base = fread( filename )
  #' Filter data and select parameters to change
  #unique(ews$part_region)
  #unique(ews$part_age_group)
  #unique(ews_base$setting)
  ews_base = subset( ews_base , ews_base$part_region       == "All" )
  ews_base = subset( ews_base , ews_base$part_gender       == "All" )
  ews_base = subset( ews_base , ews_base$part_social_group == "All" )
  ews_base = subset( ews_base , ews_base$part_income       == "All" )
  ews_base = subset( ews_base , ews_base$part_high_risk    == "All" )
  ews_base = subset( ews_base , ews_base$part_work_place   == "All" )
  ews_base = subset( ews_base , ews_base$setting           == "All" )
  
  ews_base$mid_date <- as.Date( ews_base$mid_date , format = "%d/%m/%Y") # Change format of date
  
  #"All", "0-4", "5-11", "5-17", "All-adults", "18-59", "60+", "18-29", "30-39", "40-49", "50-59", "60-69", "70+"
  ews_All = subset( ews_base , ews_base$part_age_group    == "All" )[ , c( "mid_date" , "mean" ) ]
  ews_00_04 = subset( ews_base , ews_base$part_age_group    == "0-4" )[ , c( "mid_date" , "mean" ) ] 
  ews_05_11 = subset( ews_base , ews_base$part_age_group    == "5-11" )[ , c( "mid_date" , "mean" ) ] 
  ews_05_17 = subset( ews_base , ews_base$part_age_group    == "5-17" )[ , c( "mid_date" , "mean" ) ] 
  ews_All_adults = subset( ews_base , ews_base$part_age_group    == "All-adults" )[ , c( "mid_date" , "mean" ) ] 
  ews_18_59 = subset( ews_base , ews_base$part_age_group    == "18-59" )[ , c( "mid_date" , "mean" ) ] 
  ews_60_plus = subset( ews_base , ews_base$part_age_group    == "60+" )[ , c( "mid_date" , "mean" ) ] 
  ews_18_29 = subset( ews_base , ews_base$part_age_group    == "18-29" )[ , c( "mid_date" , "mean" ) ] 
  ews_30_39 = subset( ews_base , ews_base$part_age_group    == "30-39" )[ , c( "mid_date" , "mean" ) ] 
  ews_40_49 = subset( ews_base , ews_base$part_age_group    == "40-49" )[ , c( "mid_date" , "mean" ) ] 
  ews_50_59 = subset( ews_base , ews_base$part_age_group    == "50-59" )[ , c( "mid_date" , "mean" ) ] 
  ews_60_69 = subset( ews_base , ews_base$part_age_group    == "60-69" )[ , c( "mid_date" , "mean" ) ] 
  ews_70_plus = subset( ews_base , ews_base$part_age_group    == "70+" )[ , c( "mid_date" , "mean" ) ] 
  #' Merge the new data frames
  #df_list = c(ews_All,ews_00_04,ews_05_11,ews_05_17,ews_All_adults,ews_18_29,ews_60_plus,ews_18_29,ews_30_39,ews_40_49,ews_50_59,ews_60_69,ews_70_plus)
  #ews_base = Reduce( function( x , y ) merge( x , y ,  all = TRUE ) , df_list )
  #ews_base = Reduce( dplyr::inner_join , list(ews_All,ews_00_04,ews_05_11,ews_05_17,ews_All_adults,ews_18_29,ews_60_plus,ews_18_29,ews_30_39,ews_40_49,ews_50_59,ews_60_69,ews_70_plus) )
  ews_base = merge( ews_All , ews_00_04 , by = "mid_date" , all = TRUE ) ; colnames(ews_base) = c("mid_date","All","0_4")
  ews_base = merge( ews_base , ews_05_11 , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11")
  ews_base = merge( ews_base , ews_05_17 , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11","5_17")
  ews_base = merge( ews_base , ews_All_adults , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11","5_17","All_adults")
  ews_base = merge( ews_base , ews_18_59 , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11","5_17","All_adults","18_59")
  ews_base = merge( ews_base , ews_60_plus , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11","5_17","All_adults","18_59","60_plus")
  ews_base = merge( ews_base , ews_18_29 , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11","5_17","All_adults","18_59","60_plus","18_29")
  ews_base = merge( ews_base , ews_30_39 , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11","5_17","All_adults","18_59","60_plus","18_29","30_39")
  ews_base = merge( ews_base , ews_40_49 , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11","5_17","All_adults","18_59","60_plus","18_29","30_39","40_49")
  ews_base = merge( ews_base , ews_50_59 , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11","5_17","All_adults","18_59","60_plus","18_29","30_39","40_49","50_59")
  ews_base = merge( ews_base , ews_60_69 , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11","5_17","All_adults","18_59","60_plus","18_29","30_39","40_49","50_59","60_69")
  ews_base = merge( ews_base , ews_70_plus , by = "mid_date" , all = TRUE ); colnames(ews_base) = c("mid_date","All","0_4","5_11","5_17","All_adults","18_59","60_plus","18_29","30_39","40_49","50_59","60_69","70_plus")
  
  colnames(ews_base) = c("date","CoMix_All_mean","CoMix_0_4_mean","CoMix_5_11_mean","CoMix_5_17_mean","CoMix_All_adults_mean","CoMix_18_59_mean","CoMix_60_plus_mean","CoMix_18_29_mean","CoMix_30_39_mean","CoMix_40_49_mean","CoMix_50_59_mean","CoMix_60_69_mean","CoMix_70_plus_mean")
  
  ews_base = data.frame( ews_base[ order( ews_base$date ) , ] )
  #' Because the CoMix data is weekly no leading indicator reset dates are identified
  #' Therefore need to interpolate to generate dataset of daily datapoints
  #'...
}
#' 6 - Behavioural - Google mobility
if ( lead_ind_type == 6 ){
  folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/Google Mobility'
  filename <- "Google mobility data - UK to 12 June 2022.csv"
  filename_prefix = "Google_mobility_"
  setwd( folder )
  ews_base = fread( filename )
  #' Select data to use for EWS calculation 
  ews_base = ews_base[ , c( "date" 
                            , "retail_and_recreation_percent_change_from_baseline" 
                            , "grocery_and_pharmacy_percent_change_from_baseline"
                            , "parks_percent_change_from_baseline"
                            , "transit_stations_percent_change_from_baseline"
                            , "workplaces_percent_change_from_baseline"
                            , "residential_percent_change_from_baseline"
                            ) ]
  #' Cut Google mobility data to end of April 2022 (additional data creates issues with code later)
  ews_base$date <- as.Date( ews_base$date , format = "%d/%m/%Y") # Change format of date
  ews_base = data.frame( subset( ews_base , ews_base$date <= "2022-04-30" ) )
  #' Look at reverse numbers for 'parks'
  ews_base$inv_parks_percent_change_from_baseline <- - ews_base$parks_percent_change_from_baseline
}

###############################################################################

for ( z in 1 : length( filenames ) ){
  filename = filenames[ z ]
  filename_prefix = filenames_prefix[ z ]
  setwd( folder )
  ews_base = fread( filename ) ;  ews_base[ , 1 ] = NULL ; ews_base = data.frame( ews_base )

  #' Initialise lists ready for output
  leading_indicator_dates = list()
  wave_reset_dates_lead_ind = list()
  leading_indicator_names = c()
  
  #' Loop through the various leading indicators in the leading indicator type selected
  for ( var_type in 2 : ncol( ews_base ) ){
    ################################
    #' Process EWS time series (except test case - shifted hospitalisations)
    ews = ews_base[ , c( 1 , var_type ) ] # Select date column (1) and variable column (var_type)
    colnames( ews ) <- c( "date" , "cases" ) # Rename columns
    ews = ews[ !is.na( ews$cases ) , ]# Remove days with NA values
    ews$date <- as.Date( ews$date , format = "%d/%m/%Y") # Change format of date
    ews$time <- lubridate::decimal_date( ews$date ) # Add column for decimal date
    ews$wday <- lubridate::wday( ews$date ) # Add column for day of week
    ews = ews[ Reduce( '&' , lapply( ews , is.finite ) ) , ] # Remove days with Inf values
  
    ################################
    #' Define end of a wave to reset the expanding time window for the early warning
    #' signals (EWS). Function uses method of fitting a generalised additive model (GAM)
    #' to the cases/hospitalisations (the variable we want to indicate) and 
    #' calculates the wave end/reset dates based on the first derivative of the GAM.
    #' Calculation is made on an 'add-one-in' basis, which simulates analysis in real-time
    #' The wave_reset_derivative_method() function can be found at 
    #' C:/Users/kdrake/OneDrive - Imperial College London/Documents/GitHub/Early_Warning_Signal/wave_reset_derivative_method.R
    
    #**Could do a sensitivity analysis for dates, by k and perhaps other parameters**
    #wave_reset_ix_list <- wave_reset_derivative_method( ews , dat_type) # dat_df, dat_type )
    #wave_reset_dates <- ews$date[ wave_reset_ix_list ] 
    #message( "Wave reset dates are: " , ews$date[ wave_reset_ix_list ] )
    
    #' OR
    #' Manual wave reset dates
    #' some of these dates are not in the cluster growth variance data
    #manual_wave_reset = which( dat_df$date %in% as.Date( c( "2020-08-08"
    #                                                        ,"2020-11-28"
    #                                                        ,"2021-05-06"
    #                                                        ,"2021-11-20"
    #                                                        ,"2022-02-20") ) )
    #' so use closest dates in cluster growth variance data
    #' ** still to calculate**
    #manual_wave_reset = which( dat_df$date %in% as.Date( c( "2020-08-08"
    #                                                        ,"2020-11-28"
    #                                                        ,"2021-05-06"
    #                                                        ,"2021-11-20"
    #                                                       ,"2022-02-20") ) )
    #wave_reset_ix_list <- manual_wave_reset
    #message( "Wave reset dates are: " , ews$date[ wave_reset_ix_list ] )
    #' OR use the peak of hospitalisations to reset the leading indicator statistics calculation
    #' These are rough estimates and need to be calculated
    hosp_wave_peaks = as.Date( c(  "2020-04-02"
                                  ,"2020-11-14"
                                  ,"2021-01-09"
                                  ,"2021-07-22"
                                  ,"2021-09-06"
                                  ,"2021-10-28"
                                  ,"2021-12-31"
                                  ,"2022-03-29") )
    manual_wave_reset = which( ews$date %in% hosp_wave_peaks )
    #' If not all of the hospitalisation wave peak dates are in the leading indicator
    #' data then find the closest (subsequent) dates for the resets
    if ( length( manual_wave_reset ) == length( hosp_wave_peaks ) ){
      wave_reset_ix_list <- manual_wave_reset
    } else {
      wave_reset_ix_list = c()
     for ( d in 1 : length( hosp_wave_peaks ) ){
       wave_reset_ix_list[ d ] <- which.min( abs( hosp_wave_peaks[ d ] - ews$date ) ) #'**could be improved by making the first negative value i.e. the first date after the hospitalisation wave peak**
     }
    }
    
    message( "Wave reset dates are: " , ews$date[ wave_reset_ix_list ] )
    wave_reset_dates <- ews$date[ wave_reset_ix_list ]
    
    ################################
    #' Calculate EWS indicators in expanding time window
    #' Value above cut-off level required for two consecutive periods before signal
    #' recorded. Maximum of 1 signal recorded per wave reset period 
    
    # Initialise lists
    min_window = 7 # Represents 7-time-step burn-in-period to allow values/trend to form/stabilise
    st_dev_list <- c() # replicate(min_window,NA)
    auto_cor_f_list <- c() #replicate(min_window,NA)
    s_list <- c() #replicate(min_window,NA)
    
    st_dev_signal_list <- c()
    auto_cor_f_signal_list <- c() #replicate(min_window,NA)
    s_signal_list <- c()
    
    st_dev_normalised_list <- c() #replicate(min_window,NA)
    auto_cor_f_normalised_list <- c() #replicate(min_window,NA)
    s_normalised_list <- c() #replicate(min_window,NA)
    
    ews_comp_sd_s_acf_list <- c()
    ews_comp_sd_s_list <- c()
    ews_comp_sd_acf_list <- c()
    ews_comp_s_acf_list <- c()
    
    ews_comp_sd_s_acf_signal_list <- c()
    ews_comp_sd_s_signal_list <- c()
    ews_comp_sd_acf_signal_list <- c()
    ews_comp_s_acf_signal_list <- c()
    
    # Counter for start of expanding time window (which is reset after wave)
    expanding_window_start = 1 
    
    # Loop through expanding (and resetting) time window(s)
    for (i in 1 : dim( ews )[ 1 ] ){
      
      #' Check if the current wave of leading indicator data has 'ended' and the 
      #' expanding time window needs to be reset
      if ( length( ews$date[ i - 1 ] ) == 0 ) {
        expanding_window_start = expanding_window_start
      } else {
        if ( ews$date[ i - 1 ] %in% wave_reset_dates ){
          expanding_window_start = i
        }
      }
      
      expanding_window_end <- i
      window_size <- expanding_window_end - expanding_window_start + 1
      
      # Data set within new expanding time window
      dat_window <- ews[ expanding_window_start : expanding_window_end , ]
      
      #################################  
      ###' Calculate leading indicator statistics 
      
      #' Variance (represented by standard deviation) of leading indicator
      st_dev <- sd( dat_window$cases, na.rm = TRUE ) #st_dev <- sd(dat_window$cases[!is.na(dat_window$cases)], na.rm = TRUE)
      st_dev_list <- c( st_dev_list , st_dev )
      st_dev_mean <- Reduce(  "+" 
                              , na.exclude( st_dev_list[ expanding_window_start : expanding_window_end ] ) ) / ( length( na.exclude( st_dev_list[ expanding_window_start : expanding_window_end ] ) ) )
      st_dev_st_dev <- sd( na.exclude( st_dev_list[ expanding_window_start : expanding_window_end ] ) )
      st_dev_normalised <- ( st_dev - st_dev_mean ) / st_dev_st_dev # sd(na.exclude(st_dev_list[expanding_window_start:expanding_window_end]))
      if ( length( st_dev_normalised ) == 0 ) { st_dev_normalised = NA }
      st_dev_normalised_list <- c( st_dev_normalised_list , st_dev_normalised )
      
      # Skewness
      s <- asbio::skew( ews$cases[ expanding_window_start : expanding_window_end ] ) # Alternative methods of calculating skewness: s <- prettyR::skew(ews$cases[1:expanding_window_end])$population # Returns same as asbio::skew() #s <- skewness$population(ews$cases[1:expanding_window_end])
      if( is.nan(s) || (length(s) == 0 ) ) { s = NA }
      s_list <- c( s_list , s )
      if ( length( s_list ) > 0 ){
        s_df = data.frame( s_list[ expanding_window_start : expanding_window_end ] )
        s_list_adj = s_df[ Reduce( '&' , lapply( s_df , is.finite ) ) , ]
      }
      s_mean <- Reduce( "+" , s_list_adj ) / ( length( s_list_adj ) )
      s_st_dev <- sd( s_list_adj )
      s_normalised <- ( s - s_mean ) / s_st_dev
      if ( length( s_normalised ) == 0 ) { s_normalised = NA }
      s_normalised_list <- c( s_normalised_list , s_normalised )
      
      # GAM required for calculating autocorrelation function and GAM needs a
      # minimum number of degrees of freedom (which can't be more than the number of
      # data points) to calculate a model
      if ( window_size > min_window ){
        
        #' Calculate generalised additive model (GAM) for data in expanding window 
        #' (simulates calculating indicator in real time)
        #' gam_fitting() @ C:/Users/kdrake/GitHub/Early_Warning_Signal/gam_fitting.R
        
        #' Some issues with the number of knots being greater than number of
        #' unique data points, so if there is an error then move on to next
        #' iteration in For loop.
        ews_gam_test <- try( gam_fitting(  ews[ expanding_window_start : expanding_window_end , ]
                                           , GAM_smooth_function = "cr"
                                           , deg_free_k = min( c( window_size - 5 , 50 ) ) ) #**CHECK THIS K**
                             , silent = TRUE )
        auto_cor_f_normalised = ifelse ( class( ews_gam_test[1] ) %in% 'list', 999, "Error" ) #'try-error' , "Error", NULL )
        
        if ( auto_cor_f_normalised != "Error" ){
          ews_gam <- ews_gam_test  # gam_fitting( ews[ expanding_window_start : expanding_window_end ]
          #              , GAM_smooth_function = "cr"
          #              , deg_free_k = min( c( window_size - 5 , 50 ) ) )
          # Display information for checking model quality and summary data
          #par(mfrow=c(2,2))
          #gam.check(ews_gam)
          #summary(ews_gam)
          
          # Autocorrelation at first lag (acf)
          auto_cor_f <- acf( ews_gam$residuals , plot = FALSE )$acf[ 2 ]
          auto_cor_f_list <- c( auto_cor_f_list , auto_cor_f )
          auto_cor_f_mean <- Reduce( "+" , na.exclude( auto_cor_f_list[ expanding_window_start : expanding_window_end ])) / ( length( na.exclude( auto_cor_f_list[ expanding_window_start : expanding_window_end ] ) ) )
          auto_cor_f_st_dev <- sd( na.exclude( auto_cor_f_list[ expanding_window_start : expanding_window_end ] ) )
          auto_cor_f_normalised <- ( auto_cor_f - auto_cor_f_mean ) / sd( na.exclude( auto_cor_f_list[ expanding_window_start : expanding_window_end ] ) )
        }
        
        if ( length( auto_cor_f_normalised ) == 0 ||
             auto_cor_f_normalised == "Error" ||
             auto_cor_f_normalised == 999 ||
             is.na( auto_cor_f_normalised ) ) { auto_cor_f_normalised = NA }
        
      } else { auto_cor_f_normalised = NA }
      
      auto_cor_f_normalised_list <- c( auto_cor_f_normalised_list , auto_cor_f_normalised )
      
      ##' Create composite (equal weighted) Early Warning Signals
      #' Standard deviation, skew and autocorrelation
      ews_comp_sd_s_acf <- ( st_dev_normalised + auto_cor_f_normalised + s_normalised ) / 3
      ews_comp_sd_s_acf_list <- c( ews_comp_sd_s_acf_list , ews_comp_sd_s_acf )
      #' Standard deviation, skew
      ews_comp_sd_s <- ( st_dev_normalised + s_normalised ) / 2
      ews_comp_sd_s_list <- c( ews_comp_sd_s_list , ews_comp_sd_s )
      #' Standard deviation and autocorrelation
      ews_comp_sd_acf <- ( st_dev_normalised + auto_cor_f_normalised ) / 2
      ews_comp_sd_acf_list <- c( ews_comp_sd_acf_list , ews_comp_sd_acf )
      #' Skew and autocorrelation
      ews_comp_s_acf <- ( auto_cor_f_normalised + s_normalised ) / 2
      ews_comp_s_acf_list <- c( ews_comp_s_acf_list , ews_comp_s_acf )
      
      # Record early warning signal if normalised values are greater than +2 standard deviations
      # Set signals to NA for 7 time-step (min_window) burn-in-period
      if ( window_size <= min_window ) {
        
        st_dev_signal <- NA 
        s_signal <- NA 
        auto_cor_f_signal <- NA
        
        ews_comp_sd_s_acf_signal <- NA
        ews_comp_sd_s_signal <- NA
        ews_comp_sd_acf_signal <- NA
        ews_comp_s_acf_signal <- NA
        
      } else {
        
        #' Signal recorded if two time-steps meet the 2-sigma threshold 
        #' and then not again until after wave reset (the point at which wave
        #' determined to have finished)
        #**NEED TO REDUCE REPETITION**
        st_dev_signal_list_sum <- ifelse(   is.numeric( Reduce( "+" , na.exclude( st_dev_signal_list[ expanding_window_start : expanding_window_end ] ) ) )
                                            , Reduce( "+" , na.exclude( st_dev_signal_list[ expanding_window_start : expanding_window_end ] ) )
                                            , 0 )
        auto_cor_f_signal_list_sum <- ifelse( is.numeric( Reduce( "+" , na.exclude( auto_cor_f_signal_list[ expanding_window_start : expanding_window_end ] ) ) )
                                              , Reduce( "+" , na.exclude( auto_cor_f_signal_list[ expanding_window_start : expanding_window_end ] ) )
                                              , 0 )
        s_signal_list_sum <- ifelse( is.numeric( Reduce( "+" , na.exclude( s_signal_list[ expanding_window_start : expanding_window_end ] ) ) )
                                     , Reduce( "+" , na.exclude( s_signal_list[ expanding_window_start : expanding_window_end ] ) ) 
                                     , 0 )
        
        ews_comp_sd_s_acf_signal_list_sum <- ifelse( is.numeric( Reduce( "+" , na.exclude( ews_comp_sd_s_acf_signal_list[ expanding_window_start : expanding_window_end ] ) ) )
                                                     , Reduce( "+" , na.exclude( ews_comp_sd_s_acf_signal_list[ expanding_window_start : expanding_window_end ] ) ) 
                                                     , 0 )
        
        ews_comp_sd_s_signal_list_sum <- ifelse( is.numeric( Reduce( "+" , na.exclude( ews_comp_sd_s_signal_list[ expanding_window_start : expanding_window_end ] ) ) )
                                                 , Reduce( "+" , na.exclude( ews_comp_sd_s_signal_list[ expanding_window_start : expanding_window_end ] ) ) 
                                                 , 0 )
        
        ews_comp_sd_acf_signal_list_sum <- ifelse( is.numeric( Reduce( "+" , na.exclude( ews_comp_sd_acf_signal_list[ expanding_window_start : expanding_window_end ] ) ) )
                                                   , Reduce( "+" , na.exclude( ews_comp_sd_acf_signal_list[ expanding_window_start : expanding_window_end ] ) ) 
                                                   , 0 )
        
        ews_comp_s_acf_signal_list_sum <- ifelse( is.numeric( Reduce( "+" , na.exclude( ews_comp_s_acf_signal_list[ expanding_window_start : expanding_window_end ] ) ) )
                                                  , Reduce( "+" , na.exclude( ews_comp_s_acf_signal_list[ expanding_window_start : expanding_window_end ] ) ) 
                                                  , 0 )
        
        
        if ( st_dev_signal_list_sum < 1 ) {
          st_dev_signal <- ifelse( ( st_dev_normalised_list[ i ] > 2 ) & ( st_dev_normalised_list[ i-1 ] > 2 ), 1, 0 ) 
        } else { st_dev_signal <- 0 }
        
        if (auto_cor_f_signal_list_sum < 1) {
          auto_cor_f_signal <- ifelse( ( auto_cor_f_normalised_list[ i ] > 2 ) & ( auto_cor_f_normalised_list[ i-1 ] > 2 ), 1, 0 ) 
        } else { auto_cor_f_signal <- 0 }
        
        if (s_signal_list_sum < 1) {
          s_signal <- ifelse( ( s_normalised_list[ i ] > 2) & ( s_normalised_list[ i-1 ] > 2 ), 1, 0 ) 
        } else { s_signal <- 0 }
        
        
        if (ews_comp_sd_s_acf_signal_list_sum < 1) {
          ews_comp_sd_s_acf_signal <- ifelse( ( ews_comp_sd_s_acf_list[ i ] > 2) & ( ews_comp_sd_s_acf_list[ i-1 ] > 2 ), 1, 0 ) 
        } else { ews_comp_sd_s_acf_signal <- 0 }
        
        if (ews_comp_sd_s_signal_list_sum < 1) {
          ews_comp_sd_s_signal <- ifelse( ( ews_comp_sd_s_list[ i ] > 2) & ( ews_comp_sd_s_list[ i-1 ] > 2 ), 1, 0 ) 
        } else { ews_comp_sd_s_signal <- 0 }
        
        if (ews_comp_sd_acf_signal_list_sum < 1) {
          ews_comp_sd_acf_signal <- ifelse( ( ews_comp_sd_acf_list[ i ] > 2) & ( ews_comp_sd_acf_list[ i-1 ] > 2 ), 1, 0 ) 
        } else { ews_comp_sd_acf_signal <- 0 }
        
        if (ews_comp_s_acf_signal_list_sum < 1) {
          ews_comp_s_acf_signal <- ifelse( ( ews_comp_s_acf_list[ i ] > 2) & ( ews_comp_s_acf_list[ i-1 ] > 2 ), 1, 0 ) 
        } else { ews_comp_s_acf_signal <- 0 }
        
      }
      
      #message(i," ",expanding_window_start,":",expanding_window_end," "
      #        ," Skewness= ",round(s,1)
      #        ," mean of skew= ",round(s_mean,1)
      #        ," StDev of skew= ",round(s_st_dev,1)
      #        ," normalised= ",round(s_normalised,1))#, st_dev_signal_list_sum)
      
      #' Record (binary) signals in time series lists
      st_dev_signal_list <- c( st_dev_signal_list, st_dev_signal )
      s_signal_list = c( s_signal_list, s_signal )
      auto_cor_f_signal_list = c( auto_cor_f_signal_list, auto_cor_f_signal )
      
      ews_comp_sd_s_acf_signal_list = c( ews_comp_sd_s_acf_signal_list, ews_comp_sd_s_acf_signal )
      ews_comp_sd_s_signal_list = c( ews_comp_sd_s_signal_list, ews_comp_sd_s_signal )
      ews_comp_sd_acf_signal_list = c( ews_comp_sd_acf_signal_list, ews_comp_sd_acf_signal )
      ews_comp_s_acf_signal_list = c( ews_comp_s_acf_signal_list, ews_comp_s_acf_signal )
      
      #' Define which dates meet signal threshold (via indices)
      st_dev_signal_ix <- which( st_dev_signal_list > 0 )
      s_signal_ix <- which( s_signal_list > 0 )
      auto_cor_f_signal_ix <- which( auto_cor_f_signal_list > 0 )
      
      ews_comp_sd_s_acf_signal_ix <- which( ews_comp_sd_s_acf_signal_list > 0 )
      ews_comp_sd_s_signal_ix <- which( ews_comp_sd_s_signal_list > 0 )
      ews_comp_sd_acf_signal_ix <- which( ews_comp_sd_acf_signal_list > 0 )
      ews_comp_s_acf_signal_ix <- which( ews_comp_s_acf_signal_list > 0 )
      #' Convert to indices for dat_df, which may have a different number of data points
      #' and so different dates for each index
      st_dev_signal_ix_dat <- which( dat_df$date %in% ews$date[ st_dev_signal_ix ] )
      s_signal_ix_dat <- which( dat_df$date %in% ews$date[ s_signal_ix ] )
      s_signal_ix_dat <- which( dat_df$date %in% ews$date[ s_signal_ix ] )
      auto_cor_f_signal_ix_dat <- which( dat_df$date %in% ews$date[ auto_cor_f_signal_ix ] )
      
      ews_comp_sd_s_acf_signal_ix_dat <- which( dat_df$date %in% ews$date[ ews_comp_sd_s_acf_signal_ix ] )
      ews_comp_sd_s_signal_ix_dat <- which( dat_df$date %in% ews$date[ ews_comp_sd_s_signal_ix ] )
      ews_comp_sd_acf_signal_ix_dat <- which( dat_df$date %in% ews$date[ews_comp_sd_acf_signal_ix ] )
      ews_comp_s_acf_signal_ix_dat <- which( dat_df$date %in% ews$date[ews_comp_s_acf_signal_ix ] )
    }
    ######################################
    #' Plot data
    if (i > min_window + 10 ){
      #message("chart 1: x length ",length(dat_df$date[1:expanding_window_end])," and y length ",length(dat_df$cases[1:expanding_window_end]))
      #message("chart 2: x length ",length(ews$date[1:expanding_window_end])," and y length ",length(st_dev_normalised_list))
      #message("chart 3: x length ",length(ews$date[1:expanding_window_end])," and y length ",length(s_normalised_list))
      #message("chart 4: x length ",length(ews$date[1:expanding_window_end])," and y length ",length(auto_cor_f_normalised_list))
      
      #' As the leading and lagging data don't necessarily have data on the same dates
      #' it is necessary to calculate the different plot ranges
      dat_plot_range_end <- which( dat_df$date == ews$date[ expanding_window_end ] )
      ews_plot_range_end <- expanding_window_end
      wave_reset_ix_list_dat <- which(dat_df$date %in% ews$date[wave_reset_ix_list])
      
      par(mfrow=c(3,1))
      par(xaxt="s", yaxt="s")
      par( mar = c( 5 , 5 , 2 , 2 ) ) #margins: bottom, left, top, right
      #' Hospitalisations/cases with wave reset points marked in red
      plot(dat_df$date[1:dat_plot_range_end] #expanding_window_end]
           , dat_df$cases[1:dat_plot_range_end] #expanding_window_end]
           , xlab="Date"
           , ylab=c("Covid-19",dat_type)
           , xlim=c(ews$date[1],ews$date[length(ews$date)])
           , typ = "l")
      lines(ews$date[1:expanding_window_end]
            , (ews$cases[ 1:expanding_window_end ]-min(ews$cases)) * max( dat_df$cases ) / (max( ews$cases )-min(ews$cases)) # normalising to hospitalisation data so looks ok in plot
            , xlab="Date"
            , ylab=c("Covid-19",dat_type)
            , xlim=c(ews$date[1],ews$date[length(ews$date)])
            , lty = 2
            , col= "blue")
      # test
      #lines(  dat_df$date[ expanding_window_start : expanding_window_end ]
      #      , GAM_hosp$fitted.values
      #      , col="blue"
      #      , typ="l" )
      points(dat_df$date[ wave_reset_ix_list_dat ]
             , ews$cases[ wave_reset_ix_list ]
             , col = "red"
             , cex = 3)
      for ( r in 1 : length( wave_reset_ix_list_dat ) ){
        abline( v = dat_df$date[ wave_reset_ix_list_dat[ r ] ] , lty = 2 )
      }
      legend( "topright"
              , c( dat_type , "leading indicator" , "wave reset (leading indicator", "EWS")
              , col = c( "black" , "blue" , "red", "red" )
              , cex = 1
              , lty = c(1,1,NA,NA)
              , pch = c(NA,NA,1,16))
      points( dat_df$date[ st_dev_signal_ix_dat ] , replicate( length( st_dev_signal_ix_dat ) , 10 ) , pch = 16 , col = "red" )
      points( dat_df$date[ s_signal_ix_dat ] , replicate( length( s_signal_ix_dat ) , 10 ) , pch = 16 , col = "blue" )
      points( dat_df$date[ auto_cor_f_signal_ix_dat ] , replicate( length( auto_cor_f_signal_ix_dat ) , 10 ) , pch = 16 , col = "green" )
      points( dat_df$date[ ews_comp_sd_s_acf_signal_ix_dat ] , replicate( length( ews_comp_sd_s_acf_signal_ix_dat ) , 10 ) , pch = 16 , col = "brown1" )
      points( dat_df$date[ ews_comp_sd_s_signal_ix_dat ] , replicate( length( ews_comp_sd_s_signal_ix_dat ) , 10 ) , pch = 16 , col = "cadetblue" )
      points( dat_df$date[ ews_comp_sd_acf_signal_ix_dat ] , replicate( length( ews_comp_sd_acf_signal_ix_dat ) , 10 ) , pch = 16 , col = "darkviolet" )
      points( dat_df$date[ ews_comp_s_acf_signal_ix_dat ] , replicate( length( ews_comp_s_acf_signal_ix_dat ) , 10 ) , pch = 16 , col = "orange" )
      
      
      
      
      #' Early warning signals (individual)
      plot(ews$date
           , replicate( length( ews$date ) , 0)
           , xlim = c( ews$date[1] , ews$date[ length( ews$date ) ] )
           , ylim = c(-5 , +5 )
           #    , ylim = c( min( min(   na.exclude( st_dev_normalised_list ) 
           #                          , na.exclude( auto_cor_f_normalised_list )
           #                          , na.exclude( s_normalised_list ) ) )
           #                , max( max(  na.exclude( st_dev_normalised_list )
           #                           , na.exclude( auto_cor_f_normalised_list )
           #                           , na.exclude( s_normalised_list ) ) ) )
           , xlab = "Date"
           , ylab = "Strength of Early Warning Signal (EWS)"
           , typ = "l")
      lines(  ews$date
              , replicate( length( ews$date ) , 2 ) )
      for ( r in 1:length(wave_reset_ix_list_dat)){
        abline(v = dat_df$date[wave_reset_ix_list_dat[r]], lty = 2)
      }
      legend( "bottomright"
              , c( "St.Dev." , "Skew" , "ACF")
              , col = c( "red" , "blue" , "green" )
              , cex = 1
              , lty = 1)
      #' Standard deviation
      lines( ews$date[ 1 : expanding_window_end ]
             , st_dev_normalised_list
             , col = "red")
      points( ews$date[ st_dev_signal_ix ]
              , st_dev_normalised_list[ st_dev_signal_ix ]
              , pch = 16 # filled circle marker
              , col = "red"
              , cex = 1)
      #' Skew
      lines( ews$date[ 1 : expanding_window_end ]
             , s_normalised_list
             , col = "blue")
      points(  ews$date[ s_signal_ix ]
               , s_normalised_list[ s_signal_ix ]
               , pch = 16 # filled circle marker
               , col = "blue"
               , cex = 1)
      #' ACF
      lines( ews$date[ 1:expanding_window_end ]
             , auto_cor_f_normalised_list
             , col = "green" )
      points( ews$date[ auto_cor_f_signal_ix ]
              , auto_cor_f_normalised_list[ auto_cor_f_signal_ix ]
              , pch = 16
              , col = "darkgreen"
              , cex = 1)
      
      #' Early warning signals (composite)
      plot(ews$date
           , replicate( length( ews$date ) , 0)
           , xlim = c( ews$date[1] , ews$date[ length( ews$date ) ] )
           , ylim = c(-5 , +5 )
           #    , ylim = c( min( min(   na.exclude( st_dev_normalised_list ) 
           #                          , na.exclude( auto_cor_f_normalised_list )
           #                          , na.exclude( s_normalised_list ) ) )
           #                , max( max(  na.exclude( st_dev_normalised_list )
           #                           , na.exclude( auto_cor_f_normalised_list )
           #                           , na.exclude( s_normalised_list ) ) ) )
           , xlab = "Date"
           , ylab = "Strength of Early Warning Signal (EWS)"
           , typ = "l")
      lines(  ews$date
              , replicate( length( ews$date ) , 2 ) )
      for ( r in 1:length(wave_reset_ix_list_dat)){
        abline(v = dat_df$date[wave_reset_ix_list_dat[r]], lty = 2)
      }
      legend( "bottomright"
              , c( "SD+Skew+ACF" , "SD+Skew" , "SD+ACF", "Skew+ACF")
              , col = c( "brown1" , "cadetblue" , "darkviolet", "orange" )
              , cex = 1
              , lty = 1)
      
      #' Standard deviation, Skew and ACF
      lines( ews$date[ 1 : expanding_window_end ]
             , ews_comp_sd_s_acf_list
             , col = "brown1")
      points( ews$date[ ews_comp_sd_s_acf_signal_ix ]
              , ews_comp_sd_s_acf_list[ ews_comp_sd_s_acf_signal_ix ]
              , pch = 16 # filled circle marker
              , col = "brown1"
              , cex = 1)
      #' Standard deviation & Skew
      lines( ews$date[ 1 : expanding_window_end ]
             , ews_comp_sd_s_list
             , col = "cadetblue")
      points(  ews$date[ ews_comp_sd_s_signal_ix ]
               , ews_comp_sd_s_list[ ews_comp_sd_s_signal_ix ]
               , pch = 16 # filled circle marker
               , col = "cadetblue"
               , cex = 1)
      #' Standard deviation & ACF
      lines( ews$date[ 1:expanding_window_end ]
             , ews_comp_sd_acf_list
             , col = "darkviolet" )
      points( ews$date[ ews_comp_sd_acf_signal_ix ]
              , ews_comp_sd_acf_list[ ews_comp_sd_acf_signal_ix ]
              , pch = 16
              , col = "darkviolet"
              , cex = 1)
      #' Skew & ACF
      lines( ews$date[ 1:expanding_window_end ]
             , ews_comp_s_acf_list
             , col = "orange" )
      points( ews$date[ ews_comp_s_acf_signal_ix ]
              , ews_comp_s_acf_list[ ews_comp_s_acf_signal_ix ]
              , pch = 16
              , col = "orange"
              , cex = 1)
      
      
    }
    
    #' Record EWS dates in dataframe
    l <- tibble::lst(    "Standard Deviation"         = ews$date[ st_dev_signal_ix ]
                         , "Skew"                       = ews$date[ s_signal_ix ]
                         , "Autocorrelation function"   = ews$date[ auto_cor_f_signal_ix ]
                         , "Composite: StDev+Skew+ACF"  = ews$date[ ews_comp_sd_s_acf_signal_ix ]
                         , "Composite: StDev+Skew"      = ews$date[ ews_comp_sd_s_signal_ix ]
                         , "Composite: Skew+ACF"        = ews$date[ ews_comp_s_acf_signal_ix ]
                         , "Composite: StDev+ACF"       = ews$date[ ews_comp_sd_acf_signal_ix ]
    )
  
    #' Save EWS analysis under name of particular leading indicator aspect
    new_name_df = paste( filename_prefix , colnames( ews_base )[ var_type ] , sep = "" )
    
    #' Compile EWS dates and wave reset dates under lists of dataset variables
    temp_df <- data.frame( lapply( l , `length<-`, max( lengths( l ) ) ) )
    leading_indicator_dates[[ var_type - 1 ]] = temp_df
    leading_indicator_names = c( leading_indicator_names , new_name_df )
    names( leading_indicator_dates ) <- leading_indicator_names
    wave_reset_dates_lead_ind[[ var_type - 1 ]] = wave_reset_dates
    names( wave_reset_dates_lead_ind ) <- leading_indicator_names
  
  }#' End of for loop cycling through list of various leading indicators under the particular leading indicator category chosen
  
  #' Save outputs
  if ( lead_ind_type == 2) {
    setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/analysis/p_val_filter_no/EWS calc outputs')
  } else { setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Analysis') }
  saveRDS( leading_indicator_dates , file= paste( filename_prefix , "EWS_dates.rds" , sep = "" ) )#RData
  saveRDS( leading_indicator_names , file= paste( filename_prefix , "EWS_names.rds" , sep = "" ) )#RData
  saveRDS( wave_reset_dates_lead_ind , file = paste( filename_prefix , "wave_reset_dates.rds" , sep = "" ) )#RData
}

#########################################
# Plot 
#par(mfrow=c(1,1))
#plot(ews$date,ews$cases)
#lines(ews$date[1:expanding_window_end],st_dev_normalised_list*10000+10000)

# 4 plots
par(mfrow=c(2,2))
# hospitalisations with wave reset highlighted by red points
plot( dat_df$date[ 1 : expanding_window_end ]
      , dat_df$cases[ 1 : expanding_window_end ]
      , xlab = "Date"
      , ylab = "Covid-19 hospitalisations"
      , xlim = c( ews$date[ 1 ] , ews$date[ length( ews$date ) ] ) )
#lines(dat_df$date[1:expanding_window_end]
#      , m$fitted.values[1:expanding_window_end]
#      , col="red")
#lines(ews$date[1:expanding_window_end]
#      , GAM_hosp$fitted.values[1:expanding_window_end]
#      , col="blue")
points( dat_df$date[ wave_reset_ix ]
        , dat_df$cases[ wave_reset_ix ]
        , col = "red")
# Standard deviation signal with red indicating 2sigma threshold met
plot( ews$date[ 1 : expanding_window_end ]
      , st_dev_normalised_list
      , xlim = c( ews$date[ 1 ], ews$date[ length( ews$date ) ] )
      , xlab = "date"
      , ylab = "strength of EWS signal (st.dev.)")
points( ews$date[ st_dev_signal_ix ]
        , st_dev_normalised_list[ st_dev_signal_ix ]
        , col = "red" )
lines(ews$date
      , replicate( length( ews$date ) , 2 )
      , col = "grey" )
# Skew signal with red indicating 2sigma threshold met
plot( ews$date[ 1 : expanding_window_end ]
      , s_normalised_list
      , xlim = c( ews$date[ 1 ] , ews$date[ length( ews$date ) ] )
      , xlab = "date"
      , ylab = "strength of EWS signal (skew)" )
points( ews$date[ s_signal_ix ]
        , s_normalised_list[ s_signal_ix ]
        , col = "red" )
lines( ews$date
       , replicate( length( ews$date ) , 2 ) 
       , col = "grey" )
# Auto-correlation function signal with red indicating 2sigma threshold met
plot( ews$date[ 1 : expanding_window_end ]
      , auto_cor_f_normalised_list
      , xlim = c( ews$date[ 1 ] , ews$date[ length( ews$date ) ] )
      , xlab = "date"
      , ylab = "strength of EWS signal (auto-correlation function)" )
points( ews$date[ auto_cor_f_signal_ix ]
        , auto_cor_f_normalised_list[ auto_cor_f_signal_ix ]
        , col = "red" )
lines( ews$date
       , replicate( length( ews$date ) , 2 )
       , col = "grey" )

######################################
# Plotting cases at top and signals for variety of leading indicators below
# Using Ct values (see below for normalised and transformed to viral load)
# Define data sets to plot
leading_indicator_dates = list(  hosp = ews_dates_hosp
                                 , hosp10 = ews_dates_hosp10
                                 , hosp20 = ews_dates_hosp20
                                 , hosp30 = ews_dates_hosp30
                                 , Control_Ct = ews_dates_Control_Ct
                                 , Ngene = ews_dates_Ngene
                                 , Ogene = ews_dates_Ogene
                                 , Sgene = ews_dates_Sgene
                                 , mean_Ct = ews_dates_mean_Ct
                                 , min_Ct = ews_dates_min_Ct
                                 , min_stdev = ews_dates_min_stdev
                                 , mean_stdev = ews_dates_mean_stdev
                                 , min_skew = ews_dates_min_skew
                                 , meann_skew = ews_dates_mean_skew)
leading_indicator_names = c(  "Hospitalisations"
                              , "Hospitalisations -10 days"
                              , "Hospitalisations -20 days"
                              , "Hospitalisations -30 days"
                              , "Ct control"
                              , "Ct N gene"
                              , "Ct O gene"
                              , "Ct S gene"
                              , "Ct mean of 3 genes"
                              , "Ct min of 3 genes"
                              , "Ct minimum of St.Dev. of 3 genes"
                              , "Ct mean of St.Dev. of 3 genes"
                              , "Ct minimum of skewness of 3 genes"
                              , "Ct mean of skewness of 3 genes")
plot_colour = c( "red" , "blue" , "green","brown1" , "cadetblue" , "darkviolet", "orange" )
wave_reset_dates_lead_ind = list(   hosp = wave_reset_hosp
                                    , hosp10 = wave_reset_hosp10
                                    , hosp20 = wave_reset_hosp20
                                    , hosp30 = wave_reset_hosp30
                                    , control_Ct = wave_reset_Control_Ct
                                    , Ngene = wave_reset_Ngene_Ct
                                    , Ogene = wave_reset_Ogene_Ct
                                    , Sgene = wave_reset_Sgene_Ct
                                    , mean_Ct = wave_reset_mean_Ct
                                    , min_Ct = wave_reset_min_Ct
                                    , min_stdev = wave_reset_min_stdev
                                    , mean_stdev = wave_reset_mean_stdev
                                    , min_skew = wave_reset_min_skew
                                    , mean_skew = wave_reset_mean_skew)

par( mfrow = c( 1 , 1 ) )
layout(matrix(c(1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), nrow = 21, ncol = 1, byrow = TRUE))

par( mar = c( 5 , 4 , 4 , 2 ) ) #margins: bottom, left, top, right

plot(  dat_df$date
       , dat_df$cases
       , xlab = "Date"
       , ylab = c( "Covid-19" , dat_type )
       ,las = 2
       , typ = "l" )
points(dat_df$date[ which(dat_df$date %in% wave_start_dates) ]
       , dat_df$cases[ which(dat_df$date %in% wave_start_dates) ]
       , col = "red"
       , cex = 3)
legend( "topright"
        , c( dat_type , "wave start date","St.Dev." , "Skew" , "ACF","SD+Skew+ACF" , "SD+Skew" , "SD+ACF", "Skew+ACF", "wave reset dates")
        , col = c( "black" , "red","red" , "blue" , "green","brown1" , "cadetblue" , "darkviolet", "orange" ,"black")
        , cex = 1
        , lty = c(1,NA,NA,NA,NA,NA,NA,NA,NA,2)
        , pch = c(NA,1,16,16,16,16,16,16,16,NA)
        , ncol = 3)

par(mar=c(1,4,1,2)) #margins: bottom, left, top, right

for (i in 1:length(leading_indicator_dates)){ # cycles through leading indicators
  
  temp_df = data.frame( leading_indicator_dates[ i ] )
  
  plot(temp_df[,1]
       , replicate ( length( temp_df[,1] ), 1)
       , xlim =c(min(dat_df$date),max(dat_df$date))
       , ylim = c(0,8)
       , xlab="" #Date" 
       , ylab="" # arbitrary value to separate the plots of the different statistics
       , xaxt = "n"
       , yaxt = "n"
       , col = plot_colour[1]
       , pch = 16
  )
  mtext( paste(leading_indicator_names[ i ]) , adj = 0 , padj = 0 ) #https://stackoverflow.com/questions/70155924/how-to-rotate-ylab-in-basic-plot-r
  
  temp_reset_dates = data.frame(wave_reset_dates_lead_ind[ i ])
  
  for ( r in 1:nrow(temp_reset_dates)){
    abline(v = temp_reset_dates[r,], lty = 2)
  }
  for (j in 2:7){ #7 different types of early warning signal calculated per leading indicator data type
    temp2_df = temp_df[ j ]
    points( temp2_df[,1]
            , replicate ( nrow( temp2_df ), j ) 
            , pch = 16
            , col = plot_colour[ j ]
    )
  }
}

######################################
# Plotting cases at top and signals for variety of leading indicators below
# Using normalised and transformed to viral load (see above for version with Ct values)
# Define data sets to plot
leading_indicator_dates_vl = list(#hosp = ews_dates_hosp
  #,hosp10 = ews_dates_hosp10
  #,hosp20 = ews_dates_hosp20
  #,hosp30 = ews_dates_hosp30
  #,Control_Ct = ews_dates_Control_Ct
  Ngene = ews_dates_Ngene_vl
  ,Ogene = ews_dates_Ogene_vl
  ,Sgene = ews_dates_Sgene_vl
  ,mean_vl = ews_dates_mean_vl
  ,min_vl = ews_dates_min_vl
  ,max_vl = ews_dates_max_vl
  ,min_stdev = ews_dates_min_stdev_vl
  ,mean_stdev = ews_dates_mean_stdev_vl
  ,max_stdev = ews_dates_max_stdev_vl
  ,min_skew = ews_dates_min_skew_vl
  ,mean_skew = ews_dates_mean_skew_vl
  ,max_skew = ews_dates_max_skew_vl)
leading_indicator_names_vl = c(#"Hospitalisations"
  #,"Hospitalisations -10 days"
  #,"Hospitalisations -20 days"
  #,"Hospitalisations -30 days"
  #,"Ct control"
  "N gene Ct normalised/transformed"
  ,"O gene Ct normalised/transformed"
  ,"S gene Ct normalised/transformed"
  ,"mean of 3 genes Ct normalised/transformed"
  ,"min of 3 genes Ct normalised/transformed"
  ,"max of 3 genes Ct normalised/transformed"
  ,"min of St.Dev. of 3 genes Ct normalised/transformed"
  ,"mean of St.Dev. of 3 genes Ct normalised/transformed"
  ,"max of St.Dev. of 3 genes Ct normalised/transformed"
  ,"min of skewness of 3 genes Ct normalised/transformed"
  ,"mean of skewness of 3 genes Ct normalised/transformed"
  ,"max of skewness of 3 genes Ct normalised/transformed")
plot_colour_vl = c("red" , "blue" , "green","brown1" , "cadetblue" , "darkviolet", "orange" )
wave_reset_dates_lead_ind_vl = list( #hosp = wave_reset_hosp
  #, hosp10 = wave_reset_hosp10
  #, hosp20 = wave_reset_hosp20
  #, hosp30 = wave_reset_hosp30
  #, control_Ct = wave_reset_Control_Ct
  Ngene = wave_reset_Ngene_vl
  , Ogene = wave_reset_Ogene_vl
  , Sgene = wave_reset_Sgene_vl
  , mean_vl = wave_reset_mean_vl
  , min_vl = wave_reset_min_vl
  , max_vl = wave_reset_max_vl
  , min_stdev = wave_reset_min_stdev_vl
  , mean_stdev = wave_reset_mean_stdev_vl
  , max_stdev = wave_reset_max_stdev_vl
  , min_skew = wave_reset_min_skew_vl
  , mean_skew = wave_reset_mean_skew_vl
  , max_skew = wave_reset_max_skew_vl)

par(mfrow=c(1,1))
layout(matrix(c(1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), nrow = 21, ncol = 1, byrow = TRUE))

par(mar=c(5,4,4,2)) #margins: bottom, left, top, right

plot(dat_df$date
     , dat_df$cases
     , xlab="Date"
     , ylab=c("Covid-19",dat_type)
     ,las = 2
     , typ = "l")
points(dat_df$date[ which(dat_df$date %in% wave_start_dates) ]
       , dat_df$cases[ which(dat_df$date %in% wave_start_dates) ]
       , col = "red"
       , cex = 3)
legend( "topright"
        , c( dat_type , "wave start date","St.Dev." , "Skew" , "ACF","SD+Skew+ACF" , "SD+Skew" , "SD+ACF", "Skew+ACF", "wave reset dates")
        , col = c( "black" , "red","red" , "blue" , "green","brown1" , "cadetblue" , "darkviolet", "orange" ,"black")
        , cex = 1
        , lty = c(1,NA,NA,NA,NA,NA,NA,NA,NA,2)
        , pch = c(NA,1,16,16,16,16,16,16,16,NA)
        , ncol = 3)

par(mar=c(1,4,1,2)) #margins: bottom, left, top, right

for (i in 1:length(leading_indicator_dates_vl)){ # cycles through leading indicators
  
  temp_df = data.frame( leading_indicator_dates_vl[ i ] )
  
  plot(temp_df[,1]
       , replicate ( length( temp_df[,1] ), 1)
       , xlim =c(min(dat_df$date),max(dat_df$date))
       , ylim = c(0,8)
       , xlab="" #Date" 
       , ylab="" # arbitrary value to separate the plots of the different statistics
       , xaxt = "n"
       , yaxt = "n"
       , col = plot_colour_vl[1]
       , pch = 16
  )
  mtext( paste(leading_indicator_names_vl[ i ]) , adj = 0 , padj = 0 ) #https://stackoverflow.com/questions/70155924/how-to-rotate-ylab-in-basic-plot-r
  
  temp_reset_dates = data.frame(wave_reset_dates_lead_ind_vl[ i ])
  
  for ( r in 1:nrow(temp_reset_dates)){
    abline(v = temp_reset_dates[r,], lty = 2)
  }
  for (j in 2:7){ #7 different types of early warning signal calculated per leading indicator data type
    temp2_df = temp_df[ j ]
    points( temp2_df[,1]
            , replicate ( nrow( temp2_df ), j ) 
            , pch = 16
            , col = plot_colour_vl[ j ]
    )
  }
}