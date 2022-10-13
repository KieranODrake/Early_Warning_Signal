#' Author: Kieran Drake
#' Date: May 2022
#' Purpose: Calculation of early warning signal (EWS)
#'  Based on methods described in 
#'  O’Brien, D. A., & Clements, C. F. (2021).
#'  Early warning signal reliability varies with COVID-19 waves
#'  https://doi.org/10.6084/M9.FIGSHARE.17081673.V1

################################
# install packages for fread (which is quicker than read_csv)
if(!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") 
if(!requireNamespace("Rtools", quietly = TRUE))
  install.packages("Rtools")
if(!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
library(data.table)
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

###############################
#' Functions to load
#' data_load() @ C:\Users\kdrake\OneDrive - Imperial College London\Documents\GitHub\Early_Warning_Signal
#' GAM_fitting() @ C:\Users\kdrake\OneDrive - Imperial College London\Documents\GitHub\Early_Warning_Signal
#' wave_reset_derivative_method @ C:\Users\kdrake\OneDrive - Imperial College London\Documents\GitHub\Early_Warning_Signal

###############################
#' Load case/hospitalisation data
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs'
#filename <- "data_2022-May-09 - UK Covid-19 cases by sample date.csv"
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
#dat_type <- "cases"
dat_type <- "hospitalisations"
#' data_load() @ C:\Users\kdrake\GitHub\Early_Warning_Signal
dat_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)

#' Load leading indicators

#' 1 - Test data - shifted hospitalisations
# TEMPORARILY use shifted case/hospitalisation data as leading indicator for testing
dat_df_10 <- dat_df
dat_df_20 <- dat_df
dat_df_30 <- dat_df
#dat_df_10 <- c(dat_df[-(seq(10))], rep(NA, 10))
shift <- function( x , n ){
  c( x[ -( seq( n ) ) ], rep( NA , n ) )
}
dat_df_10$cases <- shift(dat_df_10$cases, 10)
dat_df_20$cases <- shift(dat_df_20$cases, 20)
dat_df_30$cases <- shift(dat_df_30$cases, 30)

#' choose one of these as the EWS
ews = dat_df
ews <- dat_df_10
rm(dat_df_10)
ews <- dat_df_20
rm(dat_df_20)
ews <- dat_df_30
rm(dat_df_30)

#' 2 - variance of cluster logistic growth rates
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_09/Analysis'
filename <- "tfps_lgr.csv" # "tfps_lgr_wtd.csv" # "tfps_gam_lgr.csv" # "tfps_gam_lgr_wtd.csv" #"tfps_lgr_non_extant.csv" # "tfps_lgr.csv" 
setwd( folder )
ews = fread( filename )
# Select data to use for EWS calculation
#ews = ews[ , c( "date" , "ma_12_md_40" ) ] #ma = maximum age of cluster, and md = minimum descendants in cluster
ews = ews[ , c( "date" , "mina28_maxa56_md100" ) ] #mina = minimum age of cluster, maxa = maximum age of cluster and md = minimum descendants in cluster

# 3 - SARS-CoV-2 Ct values (PCR cycle threshold)
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files'
filename <- "Ct_p2_mean_df_v2.csv"
# Or 
filename <- "Ct_p2_median_df_v2.csv"
# test
setwd( folder )
ews = fread( filename )
# Select data to use for EWS calculation 
ews = ews[,c("Date","vl_max_stdev")] #O_Ct, N_Ct, S_Ct, Control_Ct, Ct_min, Ct_mean, O_Ct_norm, N_Ct_norm, S_Ct_norm, O_vl, N_vl, S_vl, vl_min, vl_mean, Ct_min_skew, Ct_min_stdev, vl_max_skew, vl_max_stdev
#dat_type = "O-gene Ct values"


#' 4 - SARS-CoV-2 positivity rates

#' 5 - Behavioural - CoMix survey
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/UK CoMix data 2022/CoMix data'
filename <- "2022-03-02_bs_means_2w_open.csv"
setwd( folder )
ews = fread( filename )
#' Filter data and select parameters to change
#unique(ews$part_region)
#unique(ews$part_age_group)
#unique(ews$setting)
ews = subset( ews , ews$part_region       == "All" )
ews = subset( ews , ews$part_gender       == "All" )
ews = subset( ews , ews$part_social_group == "All" )
ews = subset( ews , ews$part_income       == "All" )
ews = subset( ews , ews$part_high_risk    == "All" )
ews = subset( ews , ews$part_work_place   == "All" )
ews = subset( ews , ews$setting           == "All" )
ews = subset( ews , ews$part_age_group    == "18-59" ) #"All", "0-4", "5-11", "5-17", "All-adults", "18-59", "60+", "18-29", "30-39", "40-49", "50-59", "60-69", "70+"
ews = ews[ , c( "mid_date" , "mean" ) ]
ews$mid_date <- as.Date( ews$mid_date , format = "%d/%m/%Y") # Change format of date
ews = ews[ order( ews$mid_date ) , ]
#' Because the CoMix data is weekly no leading indicator reset dates are identified
#' Therefore need to interpolate to generate dataset of daily datapoints
#'...

#' 6 - Behavioural - Google mobility
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/Google Mobility'
filename <- "Google mobility data - UK to 12 June 2022.csv"
setwd( folder )
ews = fread( filename )
#' Select data to use for EWS calculation 
#"retail_and_recreation_percent_change_from_baseline"
#"grocery_and_pharmacy_percent_change_from_baseline"
#"parks_percent_change_from_baseline"
#"transit_stations_percent_change_from_baseline"
#"workplaces_percent_change_from_baseline"
#"residential_percent_change_from_baseline"
ews = ews[ , c( "date" , "retail_and_recreation_percent_change_from_baseline" ) ]
#' Cut Google mobility data to end of April 2022 (additional data creates issues with code later)
ews$date <- as.Date( ews$date , format = "%d/%m/%Y") # Change format of date
ews = subset( ews , ews$date <= "2022-04-30" )
#' Look at reverse numbers for 'parks'
ews$parks_percent_change_from_baseline <- - ews$parks_percent_change_from_baseline

################################
#' Process EWS time series (except test case - shifted hospitalisations)
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
wave_reset_ix_list <- wave_reset_derivative_method( ews , dat_type) # dat_df, dat_type )
wave_reset_dates <- ews$date[ wave_reset_ix_list ] 
message( "Wave reset dates are: " , ews$date[ wave_reset_ix_list ] )

#' OR
#' Manual wave reset dates
#' some of these dates are not in the cluster growth variance data
#manual_wave_reset = which( dat_df$date %in% as.Date( c( "2020-08-08"
#                                                        ,"2020-11-28"
#                                                        ,"2021-05-06"
#                                                        ,"2021-11-20"
#                                                        ,"2022-02-20") ) )
#' so use closest dates in cluster growth variance data
manual_wave_reset = which( ews$date %in% as.Date( c( "2020-08-14"
                                                        ,"2020-11-29"
                                                        ,"2021-05-07"
                                                        ,"2021-11-18"
                                                        ,"2022-02-21") ) )
wave_reset_ix_list <- manual_wave_reset
message( "Wave reset dates are: " , ews$date[ wave_reset_ix_list ] )
wave_reset_dates <- ews$date[ wave_reset_ix_list ]
#' OR use the peak of hospitalisations to reset the leading indicator statistics calculation
#' These are rough estimates and need to be calculated
manual_wave_reset = which( ews$date %in% as.Date( c(  "2020-04-02"
                                                     ,"2020-11-14"
                                                     ,"2021-01-09"
                                                     ,"2021-07-22"
                                                     ,"2021-09-06"
                                                     ,"2021-10-28"
                                                     ,"2021-12-31"
                                                     ,"2022-03-29") ) )
wave_reset_ix_list <- manual_wave_reset
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
    ews_gam_test <- try( gam_fitting(  ews[ expanding_window_start : expanding_window_end ]
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

#' The following is split by leading indicator type
#' 
#' 1 - Test - hospitalisations shifted
#' Convert to dataframe and fill empty cells with <NA>
#' No shift
ews_dates_hosp_shift_0 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_hosp_shift_0 = wave_reset_dates
#' Hospitalisation data shifted by 10 days
ews_dates_hosp_shift_10 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_hosp_shift_10 = wave_reset_dates
#' Hospitalisation data shifted by 20 days
ews_dates_hosp_shift_20 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_hosp_shift_20 = wave_reset_dates
#' Hospitalisation data shifted by 30 days
ews_dates_hosp_shift_30 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_hosp_shift_30 = wave_reset_dates
#' Compile into list
leading_indicator_dates = list(   hosp_shift_0 = ews_dates_hosp_shift_0
                                  , hosp_shift_10 = ews_dates_hosp_shift_10
                                  , hosp_shift_20 = ews_dates_hosp_shift_20
                                  , hosp_shift_30 = ews_dates_hosp_shift_30
)
leading_indicator_names = c(  "Hospitalisations"
                              , "Hospitalisations -10 days"
                              , "Hospitalisations -20 days"
                              , "Hospitalisations -30 days"
)

plot_colour = c( "red" , "blue" , "green","brown1" , "cadetblue" , "darkviolet", "orange" )
wave_reset_dates_lead_ind = list(   hosp_shift_0 = wave_reset_dates_hosp_shift_0
                                    , hosp_shift_10 = wave_reset_dates_hosp_shift_0
                                    , hosp_shift_20 = wave_reset_dates_hosp_shift_0
                                    , hosp_shift_30 = wave_reset_dates_hosp_shift_0
)

setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Analysis')
saveRDS( leading_indicator_dates , file="test_hosp_shift_EWS_dates.RData")
saveRDS( leading_indicator_names , file="test_hosp_shift_EWS_names.RData")
saveRDS( wave_reset_dates_lead_ind , file="test_hosp_shift_wave_reset_dates.RData")


#' 2 - Cluster logistic growth rate variance
#' Convert to dataframe and fill empty cells with <NA>
#'   #' lgrv = cluster logistic growth rate variance, ma08 = max age 8 weeks (54 days),
#'  md20 = minimum descendants in cluster, e = extant (only clusters with sequences within 14 days of tree date)
#' Need to edit these each time run with new leading indicator data
ews_dates_lgrv_ma12_md90_ne = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_lgrv_ma12_md90_ne = wave_reset_dates
#' Compile into list
leading_indicator_dates = list(   lgrv_ma08_md20_e = ews_dates_lgrv_ma08_md20_e
                                  , lgrv_ma08_md30_e = ews_dates_lgrv_ma08_md30_e
                                  , lgrv_ma08_md40_e = ews_dates_lgrv_ma08_md40_e
                                  , lgrv_ma08_md50_e = ews_dates_lgrv_ma08_md50_e
                                  , lgrv_ma08_md60_e = ews_dates_lgrv_ma08_md60_e
                                  , lgrv_ma08_md70_e = ews_dates_lgrv_ma08_md70_e
                                  , lgrv_ma08_md100_e = ews_dates_lgrv_ma08_md100_e
                                  , lgrv_ma12_md20_e = ews_dates_lgrv_ma12_md20_e
                                  , lgrv_ma12_md30_e = ews_dates_lgrv_ma12_md30_e
                                  , lgrv_ma12_md40_e = ews_dates_lgrv_ma12_md40_e
                                  , lgrv_ma12_md50_e = ews_dates_lgrv_ma12_md50_e
                                  , lgrv_ma12_md60_e = ews_dates_lgrv_ma12_md60_e
                                  , lgrv_ma12_md70_e = ews_dates_lgrv_ma12_md70_e
                                  , lgrv_ma12_md80_e = ews_dates_lgrv_ma12_md80_e
                                  , lgrv_ma12_md90_e = ews_dates_lgrv_ma12_md90_e
                                  , lgrv_ma12_md100_e = ews_dates_lgrv_ma12_md100_e
                                  , lgrv_ma08_md20_ne = ews_dates_lgrv_ma08_md20_ne
                                  , lgrv_ma08_md30_ne = ews_dates_lgrv_ma08_md30_ne
                                  , lgrv_ma08_md40_ne = ews_dates_lgrv_ma08_md40_ne
                                  , lgrv_ma08_md50_ne = ews_dates_lgrv_ma08_md50_ne
                                  , lgrv_ma08_md60_ne = ews_dates_lgrv_ma08_md60_ne
                                  , lgrv_ma08_md70_ne = ews_dates_lgrv_ma08_md70_ne
                                  , lgrv_ma08_md100_ne = ews_dates_lgrv_ma08_md100_ne
                                  , lgrv_ma12_md20_ne = ews_dates_lgrv_ma12_md20_ne
                                  , lgrv_ma12_md30_ne = ews_dates_lgrv_ma12_md30_ne
                                  , lgrv_ma12_md40_ne = ews_dates_lgrv_ma12_md40_ne
                                  , lgrv_ma12_md50_ne = ews_dates_lgrv_ma12_md50_ne
                                  , lgrv_ma12_md60_ne = ews_dates_lgrv_ma12_md60_ne
                                  , lgrv_ma12_md70_ne = ews_dates_lgrv_ma12_md70_ne
                                  , lgrv_ma12_md80_ne = ews_dates_lgrv_ma12_md80_ne
                                  , lgrv_ma12_md90_ne = ews_dates_lgrv_ma12_md90_ne
                                  , lgrv_ma12_md100_ne = ews_dates_lgrv_ma12_md100_ne
)
leading_indicator_names = c(  "Max age = 8wks, Min desc. = 20, Extant"
                              , "Max age = 8wks, Min desc. = 30, Extant"
                              , "Max age = 8wks, Min desc. = 40, Extant"
                              , "Max age = 8wks, Min desc. = 50, Extant"
                              , "Max age = 8wks, Min desc. = 60, Extant"
                              , "Max age = 8wks, Min desc. = 70, Extant"
                              , "Max age = 8wks, Min desc. = 100, Extant"
                              , "Max age = 12wks, Min desc. = 20, Extant"
                              , "Max age = 12wks, Min desc. = 30, Extant"
                              , "Max age = 12wks, Min desc. = 40, Extant"
                              , "Max age = 12wks, Min desc. = 50, Extant"
                              , "Max age = 12wks, Min desc. = 60, Extant"
                              , "Max age = 12wks, Min desc. = 70, Extant"
                              , "Max age = 12wks, Min desc. = 80, Extant"
                              , "Max age = 12wks, Min desc. = 90, Extant"
                              , "Max age = 12wks, Min desc. = 100, Non-extant"
                              , "Max age = 8wks, Min desc. = 20, Non-extant"
                              , "Max age = 8wks, Min desc. = 30, Non-extant"
                              , "Max age = 8wks, Min desc. = 40, Non-extant"
                              , "Max age = 8wks, Min desc. = 50, Non-extant"
                              , "Max age = 8wks, Min desc. = 60, Non-extant"
                              , "Max age = 8wks, Min desc. = 70, Non-extant"
                              , "Max age = 8wks, Min desc. = 100, Non-extant"
                              , "Max age = 12wks, Min desc. = 20, Non-extant"
                              , "Max age = 12wks, Min desc. = 30, Non-extant"
                              , "Max age = 12wks, Min desc. = 40, Non-extant"
                              , "Max age = 12wks, Min desc. = 50, Non-extant"
                              , "Max age = 12wks, Min desc. = 60, Non-extant"
                              , "Max age = 12wks, Min desc. = 70, Non-extant"
                              , "Max age = 12wks, Min desc. = 80, Non-extant"
                              , "Max age = 12wks, Min desc. = 90, Non-extant"
                              , "Max age = 12wks, Min desc. = 100, Non-extant"
)

plot_colour = c( "red" , "blue" , "green","brown1" , "cadetblue" , "darkviolet", "orange" )
wave_reset_dates_lead_ind = list(   lgrv_ma08_md20_e = wave_reset_dates_lgrv_ma08_md20_e
                                    , lgrv_ma08_md30_e = wave_reset_dates_lgrv_ma08_md30_e
                                    , lgrv_ma08_md40_e = wave_reset_dates_lgrv_ma08_md40_e
                                    , lgrv_ma08_md50_e = wave_reset_dates_lgrv_ma08_md50_e
                                    , lgrv_ma08_md60_e = wave_reset_dates_lgrv_ma08_md60_e
                                    , lgrv_ma08_md70_e = wave_reset_dates_lgrv_ma08_md70_e
                                    , lgrv_ma08_md100_e = wave_reset_dates_lgrv_ma08_md100_e
                                    , lgrv_ma12_md20_e = wave_reset_dates_lgrv_ma12_md20_e
                                    , lgrv_ma12_md30_e = wave_reset_dates_lgrv_ma12_md30_e
                                    , lgrv_ma12_md40_e = wave_reset_dates_lgrv_ma12_md40_e
                                    , lgrv_ma12_md50_e = wave_reset_dates_lgrv_ma12_md50_e
                                    , lgrv_ma12_md60_e = wave_reset_dates_lgrv_ma12_md60_e
                                    , lgrv_ma12_md70_e = wave_reset_dates_lgrv_ma12_md70_e
                                    , lgrv_ma12_md80_e = wave_reset_dates_lgrv_ma12_md80_e
                                    , lgrv_ma12_md90_e = wave_reset_dates_lgrv_ma12_md90_e
                                    , lgrv_ma12_md100_e = wave_reset_dates_lgrv_ma12_md100_e
                                    , lgrv_ma08_md20_ne = wave_reset_dates_lgrv_ma08_md20_ne
                                    , lgrv_ma08_md30_ne = wave_reset_dates_lgrv_ma08_md30_ne
                                    , lgrv_ma08_md40_ne = wave_reset_dates_lgrv_ma08_md40_ne
                                    , lgrv_ma08_md50_ne = wave_reset_dates_lgrv_ma08_md50_ne
                                    , lgrv_ma08_md60_ne = wave_reset_dates_lgrv_ma08_md60_ne
                                    , lgrv_ma08_md70_ne = wave_reset_dates_lgrv_ma08_md70_ne
                                    , lgrv_ma08_md100_ne = wave_reset_dates_lgrv_ma08_md100_ne
                                    , lgrv_ma12_md20_ne = wave_reset_dates_lgrv_ma12_md20_ne
                                    , lgrv_ma12_md30_ne = wave_reset_dates_lgrv_ma12_md30_ne
                                    , lgrv_ma12_md40_ne = wave_reset_dates_lgrv_ma12_md40_ne
                                    , lgrv_ma12_md50_ne = wave_reset_dates_lgrv_ma12_md50_ne
                                    , lgrv_ma12_md60_ne = wave_reset_dates_lgrv_ma12_md60_ne
                                    , lgrv_ma12_md70_ne = wave_reset_dates_lgrv_ma12_md70_ne
                                    , lgrv_ma12_md80_ne = wave_reset_dates_lgrv_ma12_md80_ne
                                    , lgrv_ma12_md90_ne = wave_reset_dates_lgrv_ma12_md90_ne
                                    , lgrv_ma12_md100_ne = wave_reset_dates_lgrv_ma12_md100_ne)

setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/Analysis')
saveRDS( leading_indicator_dates , file="lgrv_EWS_dates.RData")
saveRDS( leading_indicator_names , file="lgrv_EWS_names.RData")
saveRDS( wave_reset_dates_lead_ind , file="lgrv_wave_reset_dates.RData")


#' 3 - Ct values
#ews_dates_max_stdev_vl = data.frame( lapply( l , `length<-`, max(lengths(l))))

#' O-gene Ct value
ews_dates_p2_mean_O_Ct = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_O_Ct = wave_reset_dates
ews_dates_p2_median_O_Ct = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_O_Ct = wave_reset_dates
#' N-gene Ct value
ews_dates_p2_mean_N_Ct = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_N_Ct = wave_reset_dates
ews_dates_p2_median_N_Ct = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_N_Ct = wave_reset_dates
#' S-gene Ct value
ews_dates_p2_mean_S_Ct = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_S_Ct = wave_reset_dates
ews_dates_p2_median_S_Ct = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_S_Ct = wave_reset_dates
#' Control Ct value
ews_dates_p2_mean_Ct_control = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_Ct_control = wave_reset_dates
ews_dates_p2_median_Ct_control = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_Ct_control = wave_reset_dates
#' Daily mean of the minimum (per sample) of the three gene Ct values
ews_dates_p2_mean_Ct_min = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_Ct_min = wave_reset_dates
ews_dates_p2_median_Ct_min = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_Ct_min = wave_reset_dates
#' Daily mean of the mean (per sample) of the three gene Ct values
ews_dates_p2_mean_Ct_mean = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_Ct_mean = wave_reset_dates
ews_dates_p2_median_Ct_mean = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_Ct_mean = wave_reset_dates
#' O-gene Ct value normalised for control Ct value (gene Ct - control Ct)
ews_dates_p2_mean_O_Ct_norm = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_O_Ct_norm = wave_reset_dates
ews_dates_p2_median_O_Ct_norm = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_O_Ct_norm = wave_reset_dates
#' N-gene Ct value normalised for control Ct value (gene Ct - control Ct)
ews_dates_p2_mean_N_Ct_norm = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_N_Ct_norm = wave_reset_dates
ews_dates_p2_median_N_Ct_norm = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_N_Ct_norm = wave_reset_dates
#' S-gene Ct value normalised for control Ct value (gene Ct - control Ct)
ews_dates_p2_mean_S_Ct_norm = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_S_Ct_norm = wave_reset_dates
ews_dates_p2_median_S_Ct_norm = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_S_Ct_norm = wave_reset_dates
#' O-gene viral load proxy as per Dahdou et al (2021) = ln(2^(-(Ct - Ct_Ctrl)))
ews_dates_p2_mean_O_vl = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_O_vl = wave_reset_dates
ews_dates_p2_median_O_vl = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_O_vl = wave_reset_dates
#' N-gene viral load proxy as per Dahdou et al (2021) = ln(2^(-(Ct - Ct_Ctrl)))
ews_dates_p2_mean_N_vl = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_N_vl = wave_reset_dates
ews_dates_p2_median_N_vl = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_N_vl = wave_reset_dates
#' S-gene viral load proxy as per Dahdou et al (2021) = ln(2^(-(Ct - Ct_Ctrl)))
ews_dates_p2_mean_S_vl = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_S_vl = wave_reset_dates
ews_dates_p2_median_S_vl = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_S_vl = wave_reset_dates
#' Minimum of viral load proxy 
ews_dates_p2_mean_vl_min = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_vl_min = wave_reset_dates
ews_dates_p2_median_vl_min = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_vl_min = wave_reset_dates
#' Mean of viral load proxy 
ews_dates_p2_mean_vl_mean = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_vl_mean = wave_reset_dates
ews_dates_p2_median_vl_mean = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_vl_mean = wave_reset_dates
#' Ct min skew 
ews_dates_p2_mean_Ct_min_skew = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_Ct_min_skew = wave_reset_dates
ews_dates_p2_median_Ct_min_skew = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_Ct_min_skew = wave_reset_dates
#' Ct min stdev 
ews_dates_p2_mean_Ct_min_stdev = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_Ct_min_stdev = wave_reset_dates
ews_dates_p2_median_Ct_min_stdev = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_Ct_min_stdev = wave_reset_dates
#' viral load proxy maximum skew 
ews_dates_p2_mean_vl_max_skew = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_vl_max_skew = wave_reset_dates
ews_dates_p2_median_vl_max_skew = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_vl_max_skew = wave_reset_dates
#' viral load proxy maximum stdev 
ews_dates_p2_mean_vl_max_stdev = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_mean_vl_max_stdev = wave_reset_dates
ews_dates_p2_median_vl_max_stdev = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_p2_median_vl_max_stdev = wave_reset_dates

#' Compile into list
leading_indicator_dates = list(   p2_mean_O_Ct = ews_dates_p2_mean_O_Ct
                                , p2_mean_N_Ct = ews_dates_p2_mean_N_Ct 
                                , p2_mean_S_Ct = ews_dates_p2_mean_S_Ct
                                , p2_mean_Ct_control = ews_dates_p2_mean_Ct_control
                                , p2_mean_Ct_min = ews_dates_p2_mean_Ct_min
                                , p2_mean_Ct_mean = ews_dates_p2_mean_Ct_mean
                                , p2_mean_O_Ct_norm = ews_dates_p2_mean_O_Ct_norm
                                , p2_mean_N_Ct_norm = ews_dates_p2_mean_N_Ct_norm
                                , p2_mean_S_Ct_norm = ews_dates_p2_mean_S_Ct_norm
                                , p2_mean_O_vl = ews_dates_p2_mean_O_vl
                                , p2_mean_N_vl = ews_dates_p2_mean_N_vl
                                , p2_mean_N_vl = ews_dates_p2_mean_S_vl
                                , p2_mean_vl_min = ews_dates_p2_mean_vl_min
                                , p2_mean_vl_mean = ews_dates_p2_mean_vl_mean
                                , p2_mean_Ct_min_skew = ews_dates_p2_mean_Ct_min_skew
                                , p2_mean_Ct_min_stdev = ews_dates_p2_mean_Ct_min_stdev
                                , p2_mean_vl_max_skew = ews_dates_p2_mean_vl_max_skew
                                , p2_mean_vl_max_stdev = ews_dates_p2_mean_vl_max_stdev
                                , p2_median_O_Ct = ews_dates_p2_median_O_Ct
                                , p2_median_N_Ct = ews_dates_p2_median_N_Ct 
                                , p2_median_S_Ct = ews_dates_p2_median_S_Ct
                                , p2_median_Ct_control = ews_dates_p2_median_Ct_control
                                , p2_median_Ct_min = ews_dates_p2_median_Ct_min
                                , p2_median_Ct_mean = ews_dates_p2_median_Ct_mean
                                , p2_median_O_Ct_norm = ews_dates_p2_median_O_Ct_norm
                                , p2_median_N_Ct_norm = ews_dates_p2_median_N_Ct_norm
                                , p2_median_S_Ct_norm = ews_dates_p2_median_S_Ct_norm
                                , p2_median_O_vl = ews_dates_p2_median_O_vl
                                , p2_median_N_vl = ews_dates_p2_median_N_vl
                                , p2_median_N_vl = ews_dates_p2_median_S_vl
                                , p2_median_vl_min = ews_dates_p2_median_vl_min
                                , p2_median_vl_mean = ews_dates_p2_median_vl_mean
                                , p2_median_Ct_min_skew = ews_dates_p2_median_Ct_min_skew
                                , p2_median_Ct_min_stdev = ews_dates_p2_median_Ct_min_stdev
                                , p2_median_vl_max_skew = ews_dates_p2_median_vl_max_skew
                                , p2_median_vl_max_stdev = ews_dates_p2_median_vl_max_stdev
)
leading_indicator_names = c(  "p2_mean_O_Ct"
                            , "p2_mean_N_Ct"
                            , "p2_mean_S_Ct"
                            , "p2_mean_Ct_control"
                            , "p2_mean_Ct_min"
                            , "p2_mean_Ct_mean"
                            , "p2_mean_O_Ct_norm"
                            , "p2_mean_N_Ct_norm"
                            , "p2_mean_S_Ct_norm"
                            , "p2_mean_O_vl"
                            , "p2_mean_N_vl"
                            , "p2_mean_S_vl"
                            , "p2_mean_vl_min"
                            , "p2_mean_vl_mean"
                            , "p2_mean_Ct_min_skew"
                            , "p2_mean_Ct_min_stdev"
                            , "p2_mean_vl_max_skew"
                            , "p2_mean_vl_max_stdev"
                            , "p2_median_O_Ct"
                            , "p2_median_N_Ct"
                            , "p2_median_S_Ct"
                            , "p2_median_Ct_control"
                            , "p2_median_Ct_min"
                            , "p2_median_Ct_mean"
                            , "p2_median_O_Ct_norm"
                            , "p2_median_N_Ct_norm"
                            , "p2_median_S_Ct_norm"
                            , "p2_median_O_vl"
                            , "p2_median_N_vl"
                            , "p2_median_S_vl"
                            , "p2_median_vl_min"
                            , "p2_median_vl_mean"
                            , "p2_median_Ct_min_skew"
                            , "p2_median_Ct_min_stdev"
                            , "p2_median_vl_max_skew"
                            , "p2_median_vl_max_stdev"
)

wave_reset_dates_lead_ind = list(   p2_mean_O_Ct = wave_reset_dates_p2_mean_O_Ct
                                  , p2_mean_N_Ct = wave_reset_dates_p2_mean_N_Ct
                                  , p2_mean_S_Ct = wave_reset_dates_p2_mean_S_Ct
                                  , p2_mean_Ct_control = wave_reset_dates_p2_mean_Ct_control
                                  , p2_mean_Ct_min = wave_reset_dates_p2_mean_Ct_min
                                  , p2_mean_Ct_mean = wave_reset_dates_p2_mean_Ct_mean
                                  , p2_mean_O_Ct_norm = wave_reset_dates_p2_mean_O_Ct_norm
                                  , p2_mean_N_Ct_norm = wave_reset_dates_p2_mean_N_Ct_norm
                                  , p2_mean_S_Ct_norm = wave_reset_dates_p2_mean_S_Ct_norm
                                  , p2_mean_O_vl = wave_reset_dates_p2_mean_O_vl
                                  , p2_mean_N_vl = wave_reset_dates_p2_mean_N_vl
                                  , p2_mean_N_vl = wave_reset_dates_p2_mean_S_vl
                                  , p2_mean_vl_min = wave_reset_dates_p2_mean_vl_min
                                  , p2_mean_vl_mean = wave_reset_dates_p2_mean_vl_mean
                                  , p2_mean_Ct_min_skew = wave_reset_dates_p2_mean_Ct_min_skew
                                  , p2_mean_Ct_min_stdev = wave_reset_dates_p2_mean_Ct_min_stdev
                                  , p2_mean_vl_max_skew = wave_reset_dates_p2_mean_vl_max_skew
                                  , p2_mean_vl_max_stdev = wave_reset_dates_p2_mean_vl_max_stdev
                                  , p2_median_O_Ct = wave_reset_dates_p2_median_O_Ct
                                  , p2_median_N_Ct = wave_reset_dates_p2_median_N_Ct
                                  , p2_median_S_Ct = wave_reset_dates_p2_median_S_Ct
                                  , p2_median_Ct_control = wave_reset_dates_p2_median_Ct_control
                                  , p2_median_Ct_min = wave_reset_dates_p2_median_Ct_min
                                  , p2_median_Ct_mean = wave_reset_dates_p2_median_Ct_mean
                                  , p2_median_O_Ct_norm = wave_reset_dates_p2_median_O_Ct_norm
                                  , p2_median_N_Ct_norm = wave_reset_dates_p2_median_N_Ct_norm
                                  , p2_median_S_Ct_norm = wave_reset_dates_p2_median_S_Ct_norm
                                  , p2_median_O_vl = wave_reset_dates_p2_median_O_vl
                                  , p2_median_N_vl = wave_reset_dates_p2_median_N_vl
                                  , p2_median_N_vl = wave_reset_dates_p2_median_S_vl
                                  , p2_median_vl_min = wave_reset_dates_p2_median_vl_min
                                  , p2_median_vl_mean = wave_reset_dates_p2_median_vl_mean
                                  , p2_median_Ct_min_skew = wave_reset_dates_p2_median_Ct_min_skew
                                  , p2_median_Ct_min_stdev = wave_reset_dates_p2_median_Ct_min_stdev
                                  , p2_median_vl_max_skew = wave_reset_dates_p2_median_vl_max_skew
                                  , p2_median_vl_max_stdev = wave_reset_dates_p2_median_vl_max_stdev
)

#' Make adjustments to data frames
#leading_indicator_dates["p2_mean_Ct_min_skew"] <- ews_dates_p2_mean_Ct_min_skew
#leading_indicator_dates["p2_mean_Ct_min_st_dev"] <- ews_dates_p2_mean_Ct_min_stdev
#leading_indicator_dates["p2_mean_vl_max_skew"] <- ews_dates_p2_mean_vl_max_skew
#leading_indicator_dates["p2_mean_vl_max_st_dev"] <- ews_dates_p2_mean_vl_max_stdev#

#leading_indicator_dates["p2_median_Ct_min_skew"] <- ews_dates_p2_median_Ct_min_skew
#leading_indicator_dates["p2_median_Ct_min_st_dev"] <- ews_dates_p2_median_Ct_min_stdev
#leading_indicator_dates["p2_median_vl_max_skew"] <- ews_dates_p2_median_vl_max_skew
#leading_indicator_dates["p2_median_vl_max_st_dev"] <- ews_dates_p2_median_vl_max_stdev

#wave_reset_dates_lead_ind["p2_mean_Ct_min_skew"] <- wave_reset_dates_p2_mean_Ct_min_skew
#wave_reset_dates_lead_ind["p2_mean_Ct_min_st_dev"] <- wave_reset_dates_p2_mean_Ct_min_stdev
#wave_reset_dates_lead_ind["p2_mean_vl_max_skew"] <- wave_reset_dates_p2_mean_vl_max_skew
#wave_reset_dates_lead_ind["p2_mean_vl_max_st_dev"] <- wave_reset_dates_p2_mean_vl_max_stdev

#wave_reset_dates_lead_ind["p2_median_Ct_min_skew"] <- wave_reset_dates_p2_median_Ct_min_skew
#wave_reset_dates_lead_ind["p2_median_Ct_min_st_dev"] <- wave_reset_dates_p2_median_Ct_min_stdev
#wave_reset_dates_lead_ind["p2_median_vl_max_skew"] <- wave_reset_dates_p2_median_vl_max_skew
#wave_reset_dates_lead_ind["p2_median_vl_max_st_dev"] <- wave_reset_dates_p2_median_vl_max_stdev

#' Save data into files for loading into 'EWS_plot.R'
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Analysis')
saveRDS( leading_indicator_dates , file="Ct_EWS_dates_v2.RData")
saveRDS( leading_indicator_names , file="Ct_EWS_names_v2.RData")
saveRDS( wave_reset_dates_lead_ind , file="Ct_wave_reset_dates_v2.RData")

#' 5 - Behavioural - CoMix Survey
#' Convert to dataframe and fill empty cells with <NA>
#' all ages
ews_dates_age_all = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_all = wave_reset_dates
#' 0-4 years old
ews_dates_age_0_4 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_0_4 = wave_reset_dates
#' 5-11 years old
ews_dates_age_5_11 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_5_11 = wave_reset_dates
#' 5-17 years old
ews_dates_age_5_17 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_5_17 = wave_reset_dates
#' All-adults 
ews_dates_age_all_adults = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_all_adults = wave_reset_dates
#' 18-59 years old
ews_dates_age_18_59 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_18_59 = wave_reset_dates
#' 60+ years old
ews_dates_age_60_plus = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_60_plus = wave_reset_dates
#' 18-29 years old
ews_dates_age_18_29 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_18_29 = wave_reset_dates
#' 30-39 years old
ews_dates_age_30_39 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_30_39 = wave_reset_dates
#' 40-49 years old
ews_dates_age_40_49 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_40_49 = wave_reset_dates
#' 50-59 years old
ews_dates_age_50_59 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_50_59 = wave_reset_dates
#' 60-69 years old
ews_dates_age_60_69 = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_60_69 = wave_reset_dates
#' 70+ years old
ews_dates_age_70_plus = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_age_70_plus = wave_reset_dates


#' Compile into list
leading_indicator_dates = list(   age_all = ews_dates_age_all
                                , age_0_4 = ews_dates_age_0_4
                                , age_5_11 = ews_dates_age_5_11
                                , age_5_17 = ews_dates_age_5_17
                                , age_all_adults = ews_dates_age_all_adults
                                , age_18_59 = ews_dates_age_18_59
                                , age_60_plus = ews_dates_age_60_plus
                                , age_18_29 = ews_dates_age_18_29
                                , age_30_39 = ews_dates_age_30_39
                                , age_40_49 = ews_dates_age_40_49
                                , age_50_59 = ews_dates_age_50_59
                                , age_60_69 = ews_dates_age_60_69
                                , age_70_plus = ews_dates_age_70_plus
                                
                                
                                )
leading_indicator_names = c(  "All ages"
                              , "0-4"
                              , "5-11"
                              , "5-17"
                              , "All adults"
                              , "18-59"
                              , "60+"
                              , "18-29"
                              , "30-39"
                              , "40-49"
                              , "50-59"
                              , "60-69"
                              , "70+"
)

wave_reset_dates_lead_ind = list(   age_all = wave_reset_dates_age_all
                                    , age_0_4 = wave_reset_dates_age_0_4
                                    , age_5_11 = wave_reset_dates_age_5_11
                                    , age_5_17 = wave_reset_dates_age_5_17
                                    , age_all_adults = wave_reset_dates_age_all_adults
                                    , age_18_59 = wave_reset_dates_age_18_59
                                    , age_60_plus = wave_reset_dates_age_60_plus
                                    , age_18_29 = wave_reset_dates_age_18_29
                                    , age_30_39 = wave_reset_dates_age_30_39
                                    , age_40_49 = wave_reset_dates_age_40_49
                                    , age_50_59 = wave_reset_dates_age_50_59
                                    , age_60_69 = wave_reset_dates_age_60_69
                                    , age_60_69 = wave_reset_dates_age_60_69
                                    , age_70_plus = wave_reset_dates_age_70_plus
                                    )


setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Analysis')
saveRDS( leading_indicator_dates , file="comix_EWS_dates.RData")
saveRDS( leading_indicator_names , file="comix_EWS_names.RData")
saveRDS( wave_reset_dates_lead_ind , file="comix_wave_reset_dates.RData")


#' 6 - Google mobility
#' Convert to dataframe and fill empty cells with <NA>
#' retail_and_recreation_percent_change_from_baseline
ews_dates_retail_and_recreation = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_retail_and_recreation = wave_reset_dates
#"grocery_and_pharmacy_percent_change_from_baseline"
ews_dates_grocery_and_pharmacy = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_grocery_and_pharmacy = wave_reset_dates
#"parks_percent_change_from_baseline"
ews_dates_parks = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_parks = wave_reset_dates
# reverse "parks_percent_change_from_baseline"
ews_dates_parks_reverse = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_parks_reverse = wave_reset_dates
#"transit_stations_percent_change_from_baseline"
ews_dates_transit_stations = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_transit_stations = wave_reset_dates
#"workplaces_percent_change_from_baseline"
ews_dates_workplaces = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_workplaces = wave_reset_dates
#"residential_percent_change_from_baseline"
ews_dates_residential = data.frame( lapply( l , `length<-`, max(lengths(l))))
wave_reset_dates_residential = wave_reset_dates

#' Compile into list
leading_indicator_dates = list(   retail_and_recreation = ews_dates_retail_and_recreation
                                , grocery_and_pharmacy = ews_dates_grocery_and_pharmacy
                                , parks = ews_dates_parks
                                , parks_reverse = ews_dates_parks_reverse
                                , transit_stations = ews_dates_transit_stations
                                , workplaces = ews_dates_workplaces
                                , residential = ews_dates_residential
)
leading_indicator_names = c(  "Retail & Recreation"
                            , "Grocery & Pharmacy"
                            , "Parks"
                            , "Parks (reverse)"
                            , "Transit Stations"
                            , "Workplaces"
                            , "Residential"
)

wave_reset_dates_lead_ind = list(   retail_and_recreation = wave_reset_dates_retail_and_recreation
                                  , grocery_and_pharmacy = wave_reset_dates_grocery_and_pharmacy
                                  , parks = wave_reset_dates_parks
                                  , parks_reverse = wave_reset_dates_parks_reverse
                                  , transit_stations = wave_reset_dates_transit_stations
                                  , workplaces = wave_reset_dates_workplaces
                                  , residential = wave_reset_dates_residential
)

setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Analysis')
saveRDS( leading_indicator_dates , file="google_mobility_EWS_dates.RData")
saveRDS( leading_indicator_names , file="google_mobility_EWS_names.RData")
saveRDS( wave_reset_dates_lead_ind , file="google_mobility_wave_reset_dates.RData")

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