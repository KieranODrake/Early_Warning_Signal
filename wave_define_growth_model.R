# Adapted by KD from sc2growth.R downloaded 
# from https://gist.github.com/emvolz/3d858e03f3780b953067901121e3c7c0 on 4 May 2022

#' Estimation of case/hospitalisation growth rates with semi-parametric 
#' statistical models
#' 
#' This function estimates logistic growth rates and growth of cases 
#' from line-listed case data. 
#' Estimates are made using generalised additive models with smooth splines for
#' variation over time. 
#' Growth in cases is estimated by fitting a quasipoisson model for case numbers
#' over time with a periodic cubic spline to capture changing testing behavior
#' by day of week. 
#' 
#' @param m Generalised additive model for cases/hospitalisations generated 
#'    using mgcv::gam()
#' @param cases_df A data frame representing line-listed case data with required 
#'    columns: date (yyyy-mm-dd), number of cases/hospitalisations.
#' @param dat_type Covid-19 "cases" or "hospitalisations"
#' @param wave_bands_df Broad date ranges for UK Covid-19 waves. Used for
#'    categorising identified wave start/end dates by wave. 
#' @param GAM_smooth_function Code for smoothing function in generalised 
#'    additive model from mgcv package. See help(smooth.terms) for descriptions.
#'    Options are: "tp","ts","ds","cr","cs","sos","ps","cp","re","gp" or "so".
#' @param growth_threshold Level used to determine whether to record a start or end/reset date
#' @param deg_free_k Degrees of freedom used in GAM 
#' @return Data frame of various model information and identified wave start and end dates

# Create function
growth_method <- function( m, cases_df, dat_type, wave_bands_df, GAM_smooth_function, growth_threshold = 0, deg_free_k)
{
  #install.packages("TTR")
  library( mgcv ) 
  library( lubridate ) 
  library( glue )
  library( zoo )
  
  # Calculate smoothed log cases
  cases_df$smoothedLogCases<- as.vector(PredictMat( m$smooth[[1]], data = cases_df) %*% coef(m)[ grepl( names(coef(m)), patt = 's\\(time\\)' ) ] )

  # Organise outputs into dataframe
  growth_df = data.frame( date = cases_df$date[-1] 
  ##  , logistic_growth = diff( d$lov)
  ##  , logistic_growth_2.5 = apply( pblg, 1, function(x) quantile( x, .025 ) )
  ##  , logistic_growth_97.5 = apply( pblg, 1, function(x) quantile( x, .975 ) )
                , growth = diff( cases_df$smoothedLogCases )
                , growth_rate = TTR::ROC(cases_df$smoothedLogCases, type = "discrete")[-1]
  ##  , growth_2.5 = apply( pbg, 1, function(x) quantile( x, .025 ) )
  ##  , growth_97.5 = apply( pbg, 1, function(x) quantile( x, .975 ) )
  )
  
  # Find dates when growth goes above and below chosen threshold 
  growth_df$growth_adj = growth_df$growth - growth_threshold
  #growth_df$growth_adj = growth_df$growth_rate + threshold
  
  updn <- c(0, diff(sign(growth_df$growth_adj)))
  wave_start_ix <- which(updn > 0) #ix <- which(updn != 0)
  wave_reset_ix <- which(updn < 0)
  
  # Ensure that wave start and reset line up
  #if (length(wave_start_ix) != length(wave_reset_ix)) {
  #  ifelse(length(wave_start_ix) > length(wave_reset_ix), wave_reset_ix <- c(NA,wave_reset_ix),wave_start_ix <- c(NA,wave_start_ix))
  #}
  #if (wave_start_ix[1] > wave_reset_ix[1]){
  #  wave_start_ix <- c(NA,wave_start_ix)
  #  wave_reset_ix <- c(wave_reset_ix,NA) # assumes that the only mismatch is at the beginning of the period
  #}
  
  # Define wave start and reset dates
  wave_start_dates <- cases_df$date[wave_start_ix]
  wave_reset_dates <- cases_df$date[wave_reset_ix]
  
  # Determine whether identified wave should be retained, based on:
  # peak to trough range, peak to peak range and time resolution (time between peaks)
  # Need at least two of these rules to be breached for it to be removed as a wave
  peak_trough_cutoff = 1.379 #  Next trough must be lower than this multiplication factor
  trough_peak_cutoff = 1/peak_trough_cutoff # Next peak must be greater than this multiplication factor
  time_cutoff = 55 # days between peaks must be greater than this
  wave_start_drop_list = c()
  wave_reset_drop_list = c()
  if ( wave_reset_ix[1] < wave_start_ix[1] ){
    wave_start_ix = c(NA,wave_start_ix)
  }
  for (i in 1:length(wave_reset_ix)-1){
    peak_trough_ratio = cases_df$cases[wave_reset_ix][i-1] / cases_df$cases[wave_start_ix][i]
    trough_peak_ratio = cases_df$cases[wave_start_ix][i] / cases_df$cases[wave_reset_ix][i]
    days_peak_peak = cases_df$date[wave_reset_ix][i] - cases_df$date[wave_reset_ix][i-1]
    a = ifelse( peak_trough_ratio < peak_trough_cutoff , 1 , 0 )
    b = ifelse( trough_peak_ratio < trough_peak_cutoff , 1 , 0 )
    c = ifelse( days_peak_peak < time_cutoff , 1 , 0 )
    if( !is.na(sum(a,b,c)) ){
      if ( sum(a,b,c) == 3 ){ 
        wave_start_drop_ix = i
        wave_start_drop_list = c(wave_start_drop_list,wave_start_drop_ix)
        #wave_reset_drop_list = c(wave_reset_drop_list,wave_reset_drop)
      }
    }
  }
  
  # Trim wave start and reset dates
  if ( length(wave_start_drop_list) > 0 ){
    wave_start_adj_ix <- wave_start_ix[-wave_start_drop_list]
    wave_start_adj_dates <- cases_df$date[wave_start_adj_ix]
  }
  #wave_reset_adj_ix <- wave_reset_ix[-wave_reset_drop_list]
  #wave_reset_dates <- cases_df$date[wave_reset_ix]
  # Plot before and after adjustment for waves with low resolution
  par(mfrow=c(2,1))
  plot(cases_df$date,cases_df$cases)
  points(cases_df$date[wave_start_ix],cases_df$cases[wave_start_ix],col="red",cex=3)
  points(cases_df$date[wave_reset_ix],cases_df$cases[wave_reset_ix],col="green",cex=3)
  plot(cases_df$date,cases_df$cases)
  points(cases_df$date[wave_start_adj_ix],cases_df$cases[wave_start_adj_ix],col="red",cex=3)
  #points(cases_df$date[wave_reset_adj_ix],cases_df$cases[wave_reset_adj_ix],col="green",cex=3)
  
  
  # Return error message if more than 9 waves identified
  #if (length(wave_start_dates) > 9){
  #  message( "The maximum number of waves is 9")
  #}
  
  ## Allocate start/reset dates to particularly wave
  ## Create dataframe
  #wave_dates_banded_df = data.frame("wave_number"=c(seq(1,9,1)),"start"=c(replicate(9,NA)),"reset"=c(replicate(9,NA)))
  ## Loop through identified wave dates
  #for (i in 1:length(wave_start_dates)){
  #  # Loop through wave bands and allocate to appropriate band
  #  # If there are multiple dates within a band then the latest date will be recorded in wave_dates_banded_df
  #  for (j in 1:(length(wave_bands_df$start)-1)){
  #    
  #    wave_dates_banded_df$start[j] = ifelse(wave_start_dates[i] %in% c(wave_bands_df$start[j]:wave_bands_df$start[j+1]),
  #                                           wave_start_dates[i],
  #                                           wave_dates_banded_df$start[j])
  #  }
  #}
  
  #for (i in 1:length(wave_reset_dates)){
  #  # Loop through wave bands and allocate to appropriate band
  #  # !!!!! If there are multiple dates within a band then the latest date will be recorded in wave_dates_banded_df
  #  for (j in 1:(length(wave_bands_df$end)-1)){
  #    
  #    wave_dates_banded_df$reset[j] = ifelse(wave_reset_dates[i] %in% c(wave_bands_df$end[j]:wave_bands_df$end[j+1]),
  #                                           wave_reset_dates[i],
  #                                           wave_dates_banded_df$reset[j])
  #  }
  #}
  
  #wave_dates_banded_df$start = as.Date(wave_dates_banded_df$start, origin = "1970-01-01")
  #wave_dates_banded_df$reset = as.Date(wave_dates_banded_df$reset, origin = "1970-01-01")
  
  ##plots for testing
  #par(mfrow=c(2,1))
  #plot(cases_df$date,cases_df$cases,typ="l")
  #points(cases_df$date[wave_start_ix],cases_df$cases[wave_start_ix],col="red",cex=5)
  #points(cases_df$date[wave_reset_ix],cases_df$cases[wave_reset_ix],col="green",cex=5)
  #lines(cases_df$date,m$fitted.values,col="blue")
  
  #par(mfrow=c(1,1))
  #plot(cases_df$date,cases_df$cases,typ="l")
  #points(wave_dates_banded_df$start,replicate(length(wave_dates_banded_df$start),0),col="red",cex=5)
  #points(wave_dates_banded_df$reset,replicate(length(wave_dates_banded_df$reset),0),col="green",cex=5)
  #lines(cases_df$date,m$fitted.values,col="blue")
  
  # Capture information from gam.check()
  # In gam.check it says "Basis dimension (k) checking results. Low p-value (k-index<1)
  # may indicate that k is too low, especially if edf is close to k'."
  gam_check_capture <- capture.output(gam.check(m))
  # Sometimes k-index is in a different row and in a different position within the row
  counter =0
  for (i in 1:length(gam_check_capture)){
    line_num=c()
    # Sometimes k-index is in a different position in the row
    if (grepl( "k-index" , gam_check_capture[i] , TRUE) == TRUE){
      counter = counter + 1
      #print(c("counter:",counter))
      # counter must be equal 2 as "k-index" appears twice in the gam.check text
      if (counter == 2){
        line_num = i
        #print(c("line number",line_num))
        line_split_1 = strsplit(gam_check_capture[line_num]," ")
        line_split_2 = strsplit(gam_check_capture[line_num + 1]," ")
        if (line_num == i){
          # Sometimes k-index is a different position within the row
          ifelse ((line_split_2[[1]][7] == "") | is.na(line_split_2[[1]][7]),
                  k_index <- as.numeric(line_split_2[[1]][8],5),
                  k_index <- as.numeric(line_split_2[[1]][7],5)
          )
        }
      }
    }
  }
  k_prime <- as.numeric(strsplit(gam_check_capture[13], " ")[[1]][2],5)
  edf_value <- as.numeric(summary(m)$edf[1],5)
  
  # Define return from function
  #rv <- data.frame( "Data_type" = dat_type
  #                  , "GAM method" = summary(m)$method
  #                  , "GAM smooth function" = GAM_smooth_function
  #                  , "Degrees_of_freedom_k" = as.numeric(deg_free_k)
  #                  , "k'" = as.numeric(k_prime)
  #                  , "k-index" = as.numeric(k_index)
  #                  , "R_squared" = as.numeric(summary(m)$r.sq)
  #                  , "edf" = edf_value
  #                  , "edf/k' ratio" = edf_value/k_prime
  #                  , "GCV" = as.numeric(summary(m)$sp.criterion)
  #                  #"wave_dates" = cbind(rv_dates)
  #                  , "Wave ID method" = as.character("Growth in log(smoothed cases)")
  #                  , "Derivative_threshold" = NA
  #                  , "Derivative_significance_level" = NA
  #                  , "Growth_threshold" = growth_threshold
  #                  , "Number_of_waves" = length(wave_start_dates) #length(wave_dates_nondecimal$start),
  #                  , "Wave 1 start date" = wave_dates_banded_df$start[1]
  #                  , "Wave 2 start date" = wave_dates_banded_df$start[2]
  #                  , "Wave 3 start date" = wave_dates_banded_df$start[3]
  #                  , "Wave 4 start date" = wave_dates_banded_df$start[4]
  #                  , "Wave 5 start date" = wave_dates_banded_df$start[5]
  #                  , "Wave 6 start date" = wave_dates_banded_df$start[6]
  #                  , "Wave 7 start date" = wave_dates_banded_df$start[7]
  #                  , "Wave 8 start date" = wave_dates_banded_df$start[8]
  #                  , "Wave 9 start date" = wave_dates_banded_df$start[9]
  #                  , "Wave 1 end date" = wave_dates_banded_df$reset[1]
  #                  , "Wave 2 end date" = wave_dates_banded_df$reset[2]
  #                  , "Wave 3 end date" = wave_dates_banded_df$reset[3]
  #                  , "Wave 4 end date" = wave_dates_banded_df$reset[4]
  #                  , "Wave 5 end date" = wave_dates_banded_df$reset[5]
  #                  , "Wave 6 end date" = wave_dates_banded_df$reset[6]
  #                  , "Wave 7 end date" = wave_dates_banded_df$reset[7]
  #                  , "Wave 8 end date" = wave_dates_banded_df$reset[8]
  #                  , "Wave 9 end date" = wave_dates_banded_df$reset[9]
  #)
  
  # Define return from function
  rv <- data.frame( "Data_type" = dat_type
                    , "GAM method" = summary(m)$method
                    , "GAM smooth function" = GAM_smooth_function
                    , "Degrees_of_freedom_k" = as.numeric(deg_free_k)
                    , "k'" = as.numeric(k_prime)
                    , "k-index" = as.numeric(k_index)
                    , "R_squared" = as.numeric(summary(m)$r.sq)
                    , "edf" = edf_value
                    , "edf/k' ratio" = edf_value/k_prime
                    , "GCV" = as.numeric(summary(m)$sp.criterion)
                    #"wave_dates" = cbind(rv_dates)
                    , "Wave ID method" = as.character("Growth in log(smoothed cases)")
                    , "Derivative_threshold" = NA
                    , "Derivative_significance_level" = NA
                    , "Growth_threshold" = growth_threshold
                    , "Number_of_waves" = length(wave_start_dates) #length(wave_dates_nondecimal$start),
  )
  
  for (i in 1:100){
    j = i
    col_name = paste("wave",j,"start")
    additional_col = data.frame(placeholder_name = NA)
    names(additional_col) <- col_name
    rv = cbind(rv,additional_col[1])
  }
  for (r in 1:100){
    s = r
    col_name = paste("wave",r,"end")
    additional_col = data.frame(placeholder_name = NA)
    names(additional_col) <- col_name
    rv = cbind(rv,additional_col[1])
  }
  
  for (u in 1:length(wave_start_dates)){
    rv[u+15] = as.Date(wave_start_dates[u],origin = "1970-01-01")
  }
  
  for (v in 1:length(wave_reset_dates)){
    rv[v+15+100] = as.Date(wave_reset_dates[v],origin = "1970-01-01")
  }
  
  return(rv)
  
}