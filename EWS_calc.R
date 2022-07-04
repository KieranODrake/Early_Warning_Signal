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

###############################
#' Functions to load
#' data_load() @ C:\Users\kdrake\OneDrive - Imperial College London\Documents\GitHub\Early_Warning_Signal
#' GAM_fitting() @ C:\Users\kdrake\OneDrive - Imperial College London\Documents\GitHub\Early_Warning_Signal
#' wave_reset_derivative_method @ C:\Users\kdrake\OneDrive - Imperial College London\Documents\GitHub\Early_Warning_Signal

###############################
# Load case/hospitalisation data
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs'
#filename <- "data_2022-May-09 - UK Covid-19 cases by sample date.csv"
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
#dat_type <- "cases"
dat_type <- "hospitalisations"
# data_load() @ C:\Users\kdrake\GitHub\Early_Warning_Signal
dat_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)

# Load leading indicators
# 1 variant logistic growth rates

# 2 Ct values
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files'
filename <- "Ct_p2_mean_df.csv"
# Or 
filename <- "Ct_p2_median_df.csv"
# test
setwd( folder )
ews = fread( filename )
# Select data to use for EWS calculation 
ews = ews[,c("Date","mean_skew")] #O_Ct, N-Ct, S-Ct, Control_Ct, Ct_min, Ct_mean, min_skew, mean_skew, min_stdev, mean_stdev
#dat_type = "O-gene Ct values"
colnames(ews) <- c("date","cases") # Rename columns
ews = ews[!is.na(ews$cases),]# Remove days with NA values
ews = ews[Reduce('&', lapply(ews, is.finite)),] # Remove days with Inf values
ews$cases = 1/ews$cases
ews$date <- as.Date( ews$date , format = "%d/%m/%Y") # Change format of date
ews$time <- lubridate::decimal_date(ews$date) # Add column for decimal date
ews$wday <- lubridate::wday( ews$date ) # Add column for day of week
#ews$cases <- -ews$cases # because lower Ct value means higher viral load

#viral load test
# test
filename <- "viral_load_test.csv"
setwd( folder )
ews = fread( filename )
# Select data to use for EWS calculation 
ews = ews[,c("Date","viral load proxy")] #O_Ct, N-Ct, S-Ct, Control_Ct, Ct_min, Ct_mean, min_skew, mean_skew, min_stdev, mean_stdev
colnames(ews) <- c("date","cases") # Rename columns
ews = ews[!is.na(ews$cases),]# Remove days with NA values
ews$cases = -ews$cases
ews$date <- as.Date( ews$date , format = "%d/%m/%Y") # Change format of date
ews$time <- lubridate::decimal_date(ews$date) # Add column for decimal date
ews$wday <- lubridate::wday( ews$date ) # Add column for day of week
#ews$cases <- -ews$cases # because lower Ct value means higher viral load

# 3 positivity rates

# TEMPORARILY use shifted case/hospitalisation data as leading indicator for testing
dat_df_10 <- dat_df
dat_df_20 <- dat_df
dat_df_30 <- dat_df

# 4 Behavioural data

#dat_df_10 <- c(dat_df[-(seq(10))], rep(NA, 10))
shift <- function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}
dat_df_10$cases <- shift(dat_df_10$cases, 10)
dat_df_20$cases <- shift(dat_df_20$cases, 20)
dat_df_30$cases <- shift(dat_df_30$cases, 30)

ews = dat_df
ews <- dat_df_10
rm(dat_df_10)
ews <- dat_df_20
rm(dat_df_20)
ews <- dat_df_30
rm(dat_df_30)


################################

#' Define end of a wave to reset the expanding time window for the early warning
#' signals (EWS). Function uses method of fitting a generalised additive model (GAM)
#' to the cases/hospitalisations (the variable we want to indicate) and 
#' calculates the wave end/reset dates based on the first derivative of the GAM.
#' Calculation is made on an 'add-one-in' basis, which simulates analysis in real-time
#' The wave_reset_derivative_method() function can be found at 
#' C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/EWS calc

#**Could do a sensitivity analysis for dates, by k and perhaps other parameters**
wave_reset_ix_list <- wave_reset_derivative_method( ews, dat_type) # dat_df, dat_type )
wave_reset_dates <- ews$date[ wave_reset_ix_list ] 
message( "Wave reset dates are: ",ews$date[wave_reset_ix_list])

# Manual wave reset dates
manual_wave_reset = which(dat_df$date %in% as.Date(c("2020-08-8","2020-11-28","2021-05-06","2021-11-20","2022-02-20")))
wave_reset_ix_list <- manual_wave_reset

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
for (i in 1:dim(ews)[1]){
  
  #' Check if the current wave of leading indicator data has 'ended' and the 
  #' expanding time window needs to be reset
  if ( length( ews$date[i-1] ) == 0 ) {
    expanding_window_start = expanding_window_start
  } else {
    if ( ews$date[i-1] %in% wave_reset_dates ){
      expanding_window_start = i
    }
  }
  
  expanding_window_end <- i
  window_size <- expanding_window_end - expanding_window_start + 1

  # Data set within new expanding time window
  dat_window <- ews[expanding_window_start:expanding_window_end,]

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
  s <- asbio::skew(ews$cases[expanding_window_start:expanding_window_end]) # Alternative methods of calculating skewness: s <- prettyR::skew(ews$cases[1:expanding_window_end])$population # Returns same as asbio::skew() #s <- skewness$population(ews$cases[1:expanding_window_end])
  if( is.nan(s) || (length(s) == 0) ) { s = NA }
  s_list <- c( s_list , s )
  if (length(s_list) > 0 ){
    s_df = data.frame(s_list[expanding_window_start:expanding_window_end])
    s_list_adj = s_df[Reduce('&', lapply(s_df, is.finite)),]
  }
  s_mean <- Reduce("+", s_list_adj) / (length(s_list_adj))
  s_st_dev <- sd(s_list_adj)
  s_normalised <- (s - s_mean) / s_st_dev
  if (length(s_normalised) == 0) { s_normalised = NA }
  s_normalised_list <- c( s_normalised_list , s_normalised )
  
  # GAM required for calculating autocorrelation function and GAM needs a
  # minimum number of degrees of freedom (which can't be more than the number of
  # data points) to calculate a model
  if (window_size > min_window){
    
    #' Calculate generalised additive model (GAM) for data in expanding window (simulates calculating indicator
    #' in real time)
    #' gam_fitting() @ C:\Users\kdrake\GitHub\Early_Warning_Signal
    
    #' Some issues with the number of knots being greater than number of
    #' unique data points, so if there is an error then move on to next
    #' iteration in For loop.
    ews_gam_test <- try( gam_fitting(  ews[ expanding_window_start : expanding_window_end ]
                                     , GAM_smooth_function = "cr"
                                     , deg_free_k = min(c(window_size-5,50))) #**CHECK THIS K**
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
      auto_cor_f <- acf(ews_gam$residuals,plot=FALSE)$acf[2]
      auto_cor_f_list <- c( auto_cor_f_list , auto_cor_f )
      auto_cor_f_mean <- Reduce("+",na.exclude(auto_cor_f_list[expanding_window_start:expanding_window_end])) / (length(na.exclude(auto_cor_f_list[expanding_window_start:expanding_window_end])))
      auto_cor_f_st_dev <- sd(na.exclude(auto_cor_f_list[expanding_window_start:expanding_window_end]))
      auto_cor_f_normalised <- (auto_cor_f - auto_cor_f_mean) / sd(na.exclude(auto_cor_f_list[expanding_window_start:expanding_window_end]))
    }
    
    if (length(auto_cor_f_normalised) == 0 ||
        auto_cor_f_normalised == "Error" ||
        auto_cor_f_normalised == 999 ||
        is.na(auto_cor_f_normalised) ) { auto_cor_f_normalised = NA }
  
  } else { auto_cor_f_normalised = NA }
  
  auto_cor_f_normalised_list <- c( auto_cor_f_normalised_list , auto_cor_f_normalised )

  ##' Create composite (equal weighted) Early Warning Signals
  #' Standard deviation, skew and autocorrelation
  ews_comp_sd_s_acf <- (st_dev_normalised + auto_cor_f_normalised + s_normalised) / 3
  ews_comp_sd_s_acf_list <- c( ews_comp_sd_s_acf_list, ews_comp_sd_s_acf )
  #' Standard deviation, skew
  ews_comp_sd_s <- (st_dev_normalised + s_normalised) / 2
  ews_comp_sd_s_list <- c( ews_comp_sd_s_list, ews_comp_sd_s )
  #' Standard deviation and autocorrelation
  ews_comp_sd_acf <- (st_dev_normalised + auto_cor_f_normalised) / 2
  ews_comp_sd_acf_list <- c( ews_comp_sd_acf_list, ews_comp_sd_acf )
  #' Skew and autocorrelation
  ews_comp_s_acf <- (auto_cor_f_normalised + s_normalised) / 2
  ews_comp_s_acf_list <- c( ews_comp_s_acf_list, ews_comp_s_acf )
  
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
    
        
    if (st_dev_signal_list_sum < 1) {
      st_dev_signal <- ifelse( ( st_dev_normalised_list[ i ] > 2) & ( st_dev_normalised_list[ i-1 ] > 2 ), 1, 0 ) 
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
  st_dev_signal_ix <- which(st_dev_signal_list > 0)
  s_signal_ix <- which(s_signal_list > 0)
  auto_cor_f_signal_ix <- which(auto_cor_f_signal_list > 0)
  
  ews_comp_sd_s_acf_signal_ix <- which(ews_comp_sd_s_acf_signal_list > 0)
  ews_comp_sd_s_signal_ix <- which(ews_comp_sd_s_signal_list > 0)
  ews_comp_sd_acf_signal_ix <- which(ews_comp_sd_acf_signal_list > 0)
  ews_comp_s_acf_signal_ix <- which(ews_comp_s_acf_signal_list > 0)
  #' Convert to indices for dat_df, which may have a different number of data points
  #' and so different dates for each index
  st_dev_signal_ix_dat <- which(dat_df$date %in% ews$date[st_dev_signal_ix])
  s_signal_ix_dat <- which(dat_df$date %in% ews$date[s_signal_ix])
  s_signal_ix_dat <- which(dat_df$date %in% ews$date[s_signal_ix])
  auto_cor_f_signal_ix_dat <- which(dat_df$date %in% ews$date[auto_cor_f_signal_ix])
  
  ews_comp_sd_s_acf_signal_ix_dat <- which(dat_df$date %in% ews$date[ews_comp_sd_s_acf_signal_ix])
  ews_comp_sd_s_signal_ix_dat <- which(dat_df$date %in% ews$date[ews_comp_sd_s_signal_ix])
  ews_comp_sd_acf_signal_ix_dat <- which(dat_df$date %in% ews$date[ews_comp_sd_acf_signal_ix])
  ews_comp_s_acf_signal_ix_dat <- which(dat_df$date %in% ews$date[ews_comp_s_acf_signal_ix])
}
  ######################################
  #' Plot data
  if (i > min_window+10){
    #message("chart 1: x length ",length(dat_df$date[1:expanding_window_end])," and y length ",length(dat_df$cases[1:expanding_window_end]))
    #message("chart 2: x length ",length(ews$date[1:expanding_window_end])," and y length ",length(st_dev_normalised_list))
    #message("chart 3: x length ",length(ews$date[1:expanding_window_end])," and y length ",length(s_normalised_list))
    #message("chart 4: x length ",length(ews$date[1:expanding_window_end])," and y length ",length(auto_cor_f_normalised_list))
    
    #' As the leading and lagging data don't necessarily have data on the same dates
    #' it is necessary to calculate the different plot ranges
    dat_plot_range_end <- which(dat_df$date == ews$date[expanding_window_end])
    ews_plot_range_end <- expanding_window_end
    wave_reset_ix_list_dat <- which(dat_df$date %in% ews$date[wave_reset_ix_list])
    
    par(mfrow=c(3,1))
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
    for ( r in 1:length(wave_reset_ix_list_dat)){
      abline(v = dat_df$date[wave_reset_ix_list_dat[r]], lty = 2)
    }
    legend( "topright"
            , c( dat_type , "leading indicator" , "wave reset (leading indicator", "EWS")
            , col = c( "black" , "blue" , "red", "red" )
            , cex = 1
            , lty = c(1,1,NA,NA)
            , pch = c(NA,NA,1,16))
    points(dat_df$date[st_dev_signal_ix_dat],replicate(length(st_dev_signal_ix_dat),10),pch = 16,col="red")
    points(dat_df$date[s_signal_ix_dat],replicate(length(s_signal_ix_dat),10),pch = 16,col="blue")
    points(dat_df$date[auto_cor_f_signal_ix_dat],replicate(length(auto_cor_f_signal_ix_dat),10),pch = 16,col="green")
    points(dat_df$date[ews_comp_sd_s_acf_signal_ix_dat],replicate(length(ews_comp_sd_s_acf_signal_ix_dat),10),pch = 16,col="brown1")
    points(dat_df$date[ews_comp_sd_s_signal_ix_dat],replicate(length(ews_comp_sd_s_signal_ix_dat),10),pch = 16,col="cadetblue")
    points(dat_df$date[ews_comp_sd_acf_signal_ix_dat],replicate(length(ews_comp_sd_acf_signal_ix_dat),10),pch = 16,col="darkviolet")
    points(dat_df$date[ews_comp_s_acf_signal_ix_dat],replicate(length(ews_comp_s_acf_signal_ix_dat),10),pch = 16,col="orange")
     
    
    
    
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

# Record EWS dates in dataframe
l <- tibble::lst(  "Standard Deviation"         = ews$date[ st_dev_signal_ix ]
                   , "Skew"                       = ews$date[ s_signal_ix ]
                   , "Autocorrelation function"   = ews$date[ auto_cor_f_signal_ix ]
                   , "Composite: StDev+Skew+ACF"  = ews$date[ ews_comp_sd_s_acf_signal_ix ]
                   , "Composite: StDev+Skew"      = ews$date[ ews_comp_sd_s_signal_ix ]
                   , "Composite: Skew+ACF"        = ews$date[ ews_comp_s_acf_signal_ix ]
                   , "Composite: StDev+ACF"       = ews$date[ ews_comp_sd_acf_signal_ix ]
                )
ews_dates_mean_skew = data.frame(lapply(l, `length<-`, max(lengths(l))))

#########################################
# Plot 
#par(mfrow=c(1,1))
#plot(ews$date,ews$cases)
#lines(ews$date[1:expanding_window_end],st_dev_normalised_list*10000+10000)

# 4 plots
par(mfrow=c(2,2))
# hospitalisations with wave reset highlighted by red points
plot(dat_df$date[1:expanding_window_end]
     , dat_df$cases[1:expanding_window_end]
     , xlab="Date"
     , ylab="Covid-19 hospitalisations"
     , xlim=c(ews$date[1],ews$date[length(ews$date)]))
#lines(dat_df$date[1:expanding_window_end]
#      , m$fitted.values[1:expanding_window_end]
#      , col="red")
#lines(ews$date[1:expanding_window_end]
#      , GAM_hosp$fitted.values[1:expanding_window_end]
#      , col="blue")
points(dat_df$date[wave_reset_ix]
       , dat_df$cases[wave_reset_ix]
       , col="red")
# Standard deviation signal with red indicating 2sigma threshold met
plot(ews$date[1:expanding_window_end]
     , st_dev_normalised_list
     , xlim=c(ews$date[1], ews$date[length(ews$date)])
     , xlab = "date"
     , ylab = "strength of EWS signal (st.dev.)")
points(ews$date[st_dev_signal_ix]
       , st_dev_normalised_list[st_dev_signal_ix]
       , col="red")
lines(ews$date
      , replicate(length(ews$date),2)
      , col="grey")
# Skew signal with red indicating 2sigma threshold met
plot(ews$date[1:expanding_window_end]
     , s_normalised_list
     , xlim=c(ews$date[1], ews$date[length(ews$date)])
     ,xlab="date"
     ,ylab = "strength of EWS signal (skew)")
points(ews$date[s_signal_ix]
       , s_normalised_list[s_signal_ix]
       , col="red")
lines(ews$date
      , replicate(length(ews$date),2) 
      , col="grey")
# Auto-correlation function signal with red indicating 2sigma threshold met
plot(ews$date[1:expanding_window_end]
     , auto_cor_f_normalised_list
     , xlim=c(ews$date[1], ews$date[length(ews$date)])
     ,xlab="date"
     ,ylab = "strength of EWS signal (auto-correlation function)")
points(ews$date[auto_cor_f_signal_ix]
       , auto_cor_f_normalised_list[auto_cor_f_signal_ix]
       , col="red")
lines(ews$date
      , replicate(length(ews$date),2)
      , col="grey")

######################################
# Plotting cases at top and signals for variety of leading indicators below
# Define data sets to plot
leading_indicator_dates = list(hosp = ews_dates_hosp
                               ,hosp10 = ews_dates_hosp10
                               ,hosp20 = ews_dates_hosp20
                               ,hosp30 = ews_dates_hosp30
                               ,Control_Ct = ews_dates_Control_Ct
                               ,Ngene = ews_dates_Ngene
                               ,Ogene = ews_dates_Ogene
                               ,Sgene = ews_dates_Sgene
                               ,mean_Ct = ews_dates_mean_Ct
                               ,min_Ct = ews_dates_min_Ct
                               ,min_stdev = ews_dates_min_stdev
                               ,mean_stdev = ews_dates_mean_stdev
                               ,min_skew = ews_dates_min_skew
                               ,meann_skew = ews_dates_mean_skew)
leading_indicator_names = c("Hospitalisations"
                            ,"Hospitalisations -10 days"
                            ,"Hospitalisations -20 days"
                            ,"Hospitalisations -30 days"
                            ,"Ct control"
                            ,"Ct N gene"
                            ,"Ct O gene"
                            ,"Ct S gene"
                            ,"Ct mean of 3 genes"
                            ,"Ct min of 3 genes"
                            ,"Ct minimum of St.Dev. of 3 genes"
                            ,"Ct mean of St.Dev. of 3 genes"
                            ,"Ct minimum of skewness of 3 genes"
                            ,"Ct mean of skewness of 3 genes")
plot_colour = c("red" , "blue" , "green","brown1" , "cadetblue" , "darkviolet", "orange" )
wave_reset_dates_lead_ind = list( hosp = wave_reset_hosp
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
