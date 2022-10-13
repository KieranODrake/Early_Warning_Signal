#' Description: Function to define Covid-19 wave end point based on cases/hospitalisations.
#' Based on the approach of O'Brien & Clements (2021) 
#' https://doi.org/10.6084/m9.figshare.c.5724054
#' The approach looks at the first derivative of a generalised additive model (GAM)
#' of the cases/hospitalisations. The end of a wave is defined as the date which 
#' follows 7 days of descreasingly negative changes and 7 days of stationary
#' cases/hospitalisations.
#' Author: Kieran Drake 
#' Date: June 2022
#' 
#' @param dat_df dataframe containing: data, cases/hospitalisations, 
#'   time (date in time format), wday (weekday)
#' @return Nothing or wave reset date

wave_reset_derivative_method <- function( dat_df, dat_type ) 
{
  
  ################################
  # install required packages 
  if(!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager") 
  if(!requireNamespace("Rtools", quietly = TRUE))
    install.packages("Rtools")
  if(!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table")
  library(data.table)
  if(!requireNamespace("R.utils", quietly = TRUE))
    install.packages('R.utils')
  if(!requireNamespace("TTR", quietly = TRUE))
  install.packages("TTR")
  if(!requireNamespace("gratia", quietly = TRUE))
  install.packages("gratia") #https://gavinsimpson.github.io/gratia/
  
  ##########################
  
  # Set parameters
  min_window = 10 # Represents 10-time-step burn-in-period to allow values/trend to form/stabilise
  expanding_window_start = 1 # Counter for start of expanding time window (which is reset after wave)
  #negative_test_result = 0 # Switch for sequence of decreasing derivatives being met
  
  # Initialise lists
  wave_reset_ix_list <- c()
  
  #' Loop through dataset on an 'add-one-in' basis to simulate real-time analysis
  for (i in 1:dim(dat_df)[1]){
    
    #' Use expanding time window approach to increase size of time series, but 
    #' reset window to zero when end of wave identified
    expanding_window_end <- i
    window_size <- expanding_window_end - expanding_window_start + 1
    #' Expanding time window reset based on derivative of GAM fitted to 
    #' case/hospitalisation data.
    #' Reset test results ahead of next test
    deriv_test_window_sum = NA
    number_zeros = NA
    deriv_test_window_sum_result = 0
    stationary_deriv_test_result = 0
  
    #' Define end of a wave to reset the expanding time window for the early
    #' warning signals (EWS) Calculate generalised additive model (GAM) of
    #' cases/hospitalisations (the variable we want to indicate)
    #' gam_fitting() function @ C:\Users\kdrake\GitHub\Early_Warning_Signal
    
    #' time series must be of a minimum length for the GAM and derivative to be
    #' calculated
    if (window_size >= min_window){
      
        #' Fit GAM to expanding data window
        
        #' Some issues with the number of knots being greater than number of
        #' unique data points, so if there is an error then move on to next
        #' iteration in For loop.
        GAM_test <- try( gam_fitting(  dat_df[ expanding_window_start : expanding_window_end ]
                                     , GAM_smooth_function = "ps"
                                     , deg_free_k = min(c(window_size-5 , 30))) #**CHECK THIS K** 30 used in min test previously
                       , silent = TRUE )

        ifelse ( class( GAM_test ) %in% 'try-error' , next, "continue" )
        GAM_hosp <- gam_fitting( dat_df[ expanding_window_start : expanding_window_end ]
                                       , GAM_smooth_function = "ps"
                                       , deg_free_k = min( c( window_size-5 , 30 ) ) )#**CHECK THIS K** 30 used in min test previously
                    
        #' Calculate derivative of fitted GAM
        GAM_derivative <- gratia::derivatives(   GAM_hosp 
                                               , term = "s(time)" 
                                               , n = length( GAM_hosp$fitted.values ) 
                                               , level = 0.95 )

        #' Test for t-13 to t-7 being increasingly negative derivative
        #' and t-6 to t0 being stationary i.e. derivative = 0 within confidence interval
        if (window_size >= 14){
          test_window_start <- window_size - 14 + 1
          test_window_mid <- window_size - 7 + 1
          
          #' Ensure the derivative signs in the first half (7 days)of the test 
          #' window (14 days) are all the same and are all negative 
          deriv_test_window <- GAM_derivative$derivative[ test_window_start : test_window_mid ]
          len <- length( deriv_test_window )
          sub_len <- length( subset( deriv_test_window , deriv_test_window <=0 ) )
          
          if ( len  == sub_len ){
          
          #if ( length( unique( sign( GAM_derivative$derivative[ test_window_start : test_window_mid ] ) ) ) == 1 ){
            #if (unique( sign( GAM_derivative$derivative[ test_window_start : test_window_mid ] ) ) == -1 ){ 
              
            #' sign(diff(data)) will return a series of +1s if value is increasing,
            #' so for a series of length seven the sum should be +7
            
            deriv_test_window_sum <- Reduce( "+" , sign( diff( deriv_test_window ) ) )
              
            if ( deriv_test_window_sum == +7 ){
                
              deriv_test_window_sum_result = 1
                
              #' If derivative = zero lies within confidence interval then set 
              #' derivative to zero in order to identify stationary points in wave.
              #' If the signs for the upper and lower limits of the confidence 
              #' interval are different then zero must lie in between.
              stationary_ix <- which( sign( GAM_derivative$upper ) != sign( GAM_derivative$lower ) )
              derivative_adj <- GAM_derivative$derivative 
              derivative_adj[stationary_ix] <- 0
              stationary_test_series <- derivative_adj[ test_window_mid : window_size ]
              number_zeros <- length( subset( stationary_test_series, stationary_test_series == 0 ) )
                
              #' Require period of 7 days with stationary change in cases (derivative = 0 for 7 days)
                #if ( length( unique( stationary_test_series ) ) == 1 ){
                #  if ( unique( stationary_test_series ) == 0 ) {
              if ( number_zeros == 7 ){  
                stationary_deriv_test_result = 1
                wave_reset_ix <- i
                wave_reset_ix_list <- c(wave_reset_ix_list, wave_reset_ix)
                expanding_window_start <- i
              }
            }
          }
          
          #' Printing messages for testing purposes
          #message( i, " "
          #        , ifelse(!exists("deriv_test_window_sum"),"temp",deriv_test_window_sum), " "
          #        , ifelse(!exists("number_zeros"),"temp",number_zeros), " "
          #        , round(GAM_derivative$derivative[ test_window_start : test_window_mid ]), " "
          #        , ifelse(!exists("stationary_test_series"),"temp",round(stationary_test_series)), " "
          #        )
        }
    }
  }
  
  #' Plot dataset with wave reset points marked with red circles
  #' If plotting on 'add-one-in' basis is required then copy and paste in previous loop
  par(mfrow=c(1,1))
  plot(  dat_df$time[1:i]
         , dat_df$cases[1:i]
         , typ="l"
         , xlab = "date"
         , ylab = c("Covid-19 ",dat_type))
  #lines(  dat_df$time[expanding_window_start : expanding_window_end]
  #        , dat_df$cases[expanding_window_start : expanding_window_end]
  #        , typ="l"
  #        , col="blue")
  points(  dat_df$time[wave_reset_ix_list]
           , dat_df$cases[wave_reset_ix_list] 
           , col="red"
           , cex=3)
  
  return(wave_reset_ix_list)
}