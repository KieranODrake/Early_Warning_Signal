#' Plot of days between early warning signal (EWS) dates for different leading 
#' indicators and the 'start date' of the Covid-19 hospitalisation wave.
#' EWS dates calculated using 
#' EWS_calc.R
#' Wave start dates calculated using
#' wave_define_growth_model.R and wave - define - hosp - growth - overlay.R
#' All in C:/Users/kdrake/OneDrive - Imperial College London/Documents/GitHub/Early_Warning_Signal/
#' 

#' 1 - Load data
#' Wave start dates (average where a wave has more than one start date as 
#' produced from the 12 models incorporating 'low resolution wave filter' - see
#' 21 June 2022 update slides)
wave_start_dates <- c( as.Date("")           #' 1 Wuhan
                      ,as.Date("2020-08-19") #' 2 B.1.177
                      ,as.Date("2020-11-29") #' 3 Alpha
                      ,as.Date("2021-05-11") #' 4 Delta
                      ,as.Date("2021-08-03") #' 5 Delta
                      ,as.Date("2021-09-27") #' 6 Delta
                      ,as.Date("2021-11-26") #' 7 Omicron
                      ,as.Date("2022-02-21") #' 8 Omicron
                      ,as.Date("")           #' 9 Omicron
                      )
#' Load EWS dates (choose type of leading indicator to work with)
#' 1.1 - Test - shifted hospitalisation data
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Analysis')
ews_dates <- readRDS( "test_hosp_shift_EWS_dates.rds" )
ews_names <- readRDS( "test_hosp_shift_EWS_names.rds" )
wave_reset_dates <- readRDS( "test_hosp_shift_wave_reset_dates.rds")
#' 1.2.1 - Cluster logistic growth rate variance
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_09/Analysis')
#"vlgr_mdperc" #"vlgr_" "v_gam_lgr_" "vlgr_wtd" "lgr_max"
filename_prefix = "vlgr_mdperc_"
ews_dates        <- readRDS( paste( filename_prefix , "EWS_dates.rds" , sep = "" ) ) # EWS from cluster logistic growth rate leading indicator calculated using EWS_calc.R
ews_names        <- readRDS( paste( filename_prefix , "EWS_names.rds" , sep = "" ) )
wave_reset_dates <- readRDS( paste( filename_prefix , "wave_reset_dates.rds" , sep = "" ) ) # Cluster logistic growth rate leading indicator wave reset dates calculated using wave_reset_derivative_method.R in EWS_calc.R
#' 1.2.2 - Composite of cluster logistic growth rates leading indicators
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_09/Analysis')
ews_dates        <- readRDS( "lead_ind_comp_EWS_dates.RData" ) # EWS from cluster logistic growth rate leading indicator calculated using EWS_calc.R
ews_names        <- readRDS( "lead_ind_comp_EWS_names.RData" )
wave_reset_dates <- readRDS( "lead_ind_comp_wave_reset_dates.RData" ) # Cluster logistic growth rate leading indicator wave reset dates calculated using wave_reset_derivative_method.R in EWS_calc.R
#' 1.3 - PCR cycle threshold (Ct) values EWS
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Analysis')
ews_dates <- readRDS( "Ct_p2_mean_EWS_dates.rds")
ews_names <- readRDS( "Ct_p2_mean_EWS_names.rds")
wave_reset_dates <- readRDS( "Ct_p2_mean_wave_reset_dates.rds")
#' Or
ews_dates <- readRDS( "Ct_p2_median_EWS_dates.rds")
ews_names <- readRDS( "Ct_p2_median_EWS_names.rds")
wave_reset_dates <- readRDS( "Ct_p2_median_wave_reset_dates.rds")
#' 1.4 - PCR positivity rate in England
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_09/Analysis')
#"vlgr_mdperc" #"vlgr_" "v_gam_lgr_" "vlgr_wtd" "lgr_max"
filename_prefix = "PCR_positivity_rate_"
ews_dates        <- readRDS( paste( filename_prefix , "EWS_dates.rds" , sep = "" ) ) # EWS from cluster logistic growth rate leading indicator calculated using EWS_calc.R
ews_names        <- readRDS( paste( filename_prefix , "EWS_names.rds" , sep = "" ) )
wave_reset_dates <- readRDS( paste( filename_prefix , "wave_reset_dates.rds" , sep = "" ) ) # Cluster logistic growth rate leading indicator wave reset dates calculated using wave_reset_derivative_method.R in EWS_calc.R
#' 1.5 - Behavioural - CoMix Survey 
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Analysis')
filename_prefix = "CoMix_"
ews_dates = readRDS( "CoMix_EWS_dates.rds" )
ews_names = readRDS( "CoMix_EWS_names.rds" )
wave_reset_dates = readRDS( "CoMix_wave_reset_dates.rds" )
#' 1.6 - Google mobility
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Analysis')
ews_dates <- readRDS( "Google_mobility_EWS_dates.rds" )
ews_names <- readRDS( "Google_mobility_EWS_names.rds" )
wave_reset_dates <- readRDS( "Google_mobility_wave_reset_dates.rds" )

#' Manually determined start and end dates for waves to enable EWS dates to be assigned to a wave
# Manually define wave date bands 
# (semi-arbitrary earliest possible start date for wave i.e. just after previous wave peak) 
#'**1 - not overlapping so miss some EWS dates and not assigned to a wave**
wave_bands_start <- c(as.Date("2020-01-01"), # 1  Wuhan
                      as.Date("2020-04-23"), # 2  Alpha
                      as.Date("2020-11-14"), # 3  Alpha
                      as.Date("2021-01-06"), # 4  Delta
                      as.Date("2021-07-18"), # 5  Delta
                      as.Date("2021-09-05"), # 6  Delta
                      as.Date("2021-10-22"), # 7  Omicron
                      # as.Date("2022-01-05"), # 8  Omicron Not an obvious separate peak in hospitalisation data
                      as.Date("2022-01-27"), # 9  Omicron
                      as.Date("2022-03-25")  # 10 Omicron
)
# (semi-arbitrary earliest possible end date for wave i.e. mid-point of next wave peak)
wave_bands_end <- c(as.Date("2020-03-27"), # 1  Wuhan
                    as.Date("2020-10-10"), # 2  Alpha
                    as.Date("2020-12-18"), # 3  Alpha
                    as.Date("2021-06-29"), # 4  Delta
                    as.Date("2021-08-21"), # 5  Delta
                    as.Date("2021-10-11"), # 6  Delta
                    as.Date("2021-12-22"), # 7  Omicron
                    #as.Date("2022-01-18"), # 8  Omicron Not an obvious separate peak in hospitalisation data
                    as.Date("2022-03-07"), # 9  Omicron
                    as.Date("2022-03-29")  # 10 Omicron
)

#' Alternative method for defining wave bands as some EWS dates not falling witin bands
#' Manually chosen date roughly half way between peak and trough on the wave decline
#'**2**  
wave_bands_start <- c(as.Date("2020-01-01"), # 1  Wuhan
                      as.Date("2020-05-01"), # 2  Alpha
                      as.Date("2020-11-23"), # 3  Alpha
                      as.Date("2021-02-01"), # 4  Delta
                      as.Date("2021-08-01"), # 5  Delta
                      as.Date("2021-09-18"), # 6  Delta
                      as.Date("2021-11-11"), # 7  Omicron
                      # as.Date("2022-01-05"), # 8  Omicron Not an obvious separate peak in hospitalisation data
                      as.Date("2022-01-29"), # 9  Omicron
                      as.Date("2022-04-23")  # 10 Omicron
)
wave_bands_end <- c(as.Date("2020-04-30"), # 1  Wuhan
                    as.Date("2020-11-22"), # 2  Alpha
                    as.Date("2021-01-31"), # 3  Alpha
                    as.Date("2021-07-30"), # 4  Delta
                    as.Date("2021-09-17"), # 5  Delta
                    as.Date("2021-11-10"), # 6  Delta
                    as.Date("2022-01-28"), # 7  Omicron
                    #as.Date("2022-01-18"), # 8  Omicron Not an obvious separate peak in hospitalisation data
                    as.Date("2022-04-22"), # 9  Omicron
                    as.Date("2022-07-30")  # 10 Omicron
)

wave_bands <- data.frame(  "wave" = c(1,2,3,4,5,6,7,8,9)
                         , "band_start" = wave_bands_start
                         , "band_end" = wave_bands_end
                         , "wave_start" = wave_start_dates
                         )
rm( wave_bands_start , wave_bands_end, wave_start_dates)

#' Lead/lag time (in days) = Wave start dates - EWS dates. 
#' Assumption must be made as to which wave the EWS date is for.

#########################################
#' Convert ews_dates to an array. 
#' 
#' 1 - Test - hospilisations shifted
#' Dim = li_n leading indicators, 9 waves and 7 leading indicator statistics. 
#' Data in the array is the EWS date 
li_n = 4 #' Number of leading indicators analysing
w_n = 9 #' Number of waves
s_n = 7 #' Number of statistics calculated for leading indicators

#' 2 - Cluster logistic growth rate variance
#' Dim = 32 leading indicators, 9 waves and 7 
#' leading indicator statistics. Data in the array is the EWS date 
li_n = 18 #' Number of leading indicators analysing
w_n = 9 #' Number of waves
s_n = 7 #' Number of statistics calculated for leading indicators

#' 3 - Ct values
#' Dim = li_n leading indicators, 9 waves and 7 leading indicator statistics. 
#' Data in the array is the EWS date 
li_n = 19 #' Number of leading indicators analysing
w_n = 9 #' Number of waves
s_n = 7 #' Number of statistics calculated for leading indicators

#' 4 - PCR positivity rates pillar 2 in England
#' Dim = li_n leading indicators, 9 waves and 7 leading indicator statistics. 
#' Data in the array is the EWS date 
li_n = 2 #' Number of leading indicators analysing
w_n = 9 #' Number of waves
s_n = 7 #' Number of statistics calculated for leading indicators

#' 5 - Behavioural - CoMix Survey
#' Dim = li_n leading indicators, 9 waves and 7 leading indicator statistics. 
#' Data in the array is the EWS date 
li_n = 13 #' Number of leading indicators analysing
w_n = 9 #' Number of waves
s_n = 7 #' Number of statistics calculated for leading indicators

#' 6 - Behavioural - Google mobility
#' Dim = li_n leading indicators, 9 waves and 7 leading indicator statistics. 
#' Data in the array is the EWS date 
li_n = 7 #' Number of leading indicators analysing
w_n = 9 #' Number of waves
s_n = 7 #' Number of statistics calculated for leading indicators

#' Below can be used for all leading indicators to transform EWS dates from lists to an array ( number of variables , n waves , n statistics )
ews_dates_array <- array( data = NA, dim = c( li_n , w_n , s_n ))
ews_diff_array <- array( data = NA, dim = c( li_n , w_n , s_n ))

multi_date_log_list <- list()
multi_date_log = c()
number_counter = 0
number_counter_2 = 0

for ( lead_ind_i in 1 : length( ews_dates ) ) {
  
  for ( wave_j in 1 : w_n ){
    band_start <- wave_bands$band_start[ wave_j ]
    band_end   <- wave_bands$band_end[ wave_j ]
    
    for ( stat_k in 1 : 7 ) {
      
      for ( date_m in 1 : length( ews_dates[[ lead_ind_i ]][ , stat_k ] ) ){
        d = ews_dates[[ lead_ind_i ]][ date_m , stat_k ]
        
        if ( !is.na( d ) ){
          if ( ( d >= band_start ) && ( d <= band_end ) ){ #' Check if EWS date is within the wave band
            number_counter = number_counter + 1
            message( number_counter )
            
            if ( is.na( ews_dates_array[ lead_ind_i , wave_j , stat_k ] ) ){
              number_counter_2 = number_counter_2 + 1
              message(  number_counter_2 )
              if (number_counter_2 == 21){
                message( paste( lead_ind_i , wave_j , stat_k ) )
              }
              ews_dates_array[ lead_ind_i , wave_j , stat_k ] = as.Date( d , origin = "1970-01-01" ) 
              ews_diff_array[ lead_ind_i , wave_j , stat_k ] = wave_bands$wave_start[ wave_j ] - d
              message( paste( number_counter_2 , sum( !is.na( ews_dates_array[  , , ] )) ) )
            } else {
              multi_date_log <- c( lead_ind_i , wave_j , stat_k )
              multi_date_log_list <- rbind( multi_date_log_list , multi_date_log )
            }
          }
        }
      }
    }
  }
}

#' Test that array conversion has worked

#' 1 - Test leading indicator - hospitalisations shifted
#'#' Dim = 32 leading indicators, 9 waves and 7 leading indicator statistics. 
#' Data in the array is the EWS date 
as.Date(ews_dates_array[ 1 , 2 , 1 ],origin="1970-01-01")
ews_dates_array[ 1 , 2 , 1 ] == ews_dates[[1]][ 1 , 1]
as.Date(ews_dates_array[ 1 , 8 , 1 ],origin="1970-01-01")
ews_dates_array[ 1 , 8 , 1 ] == ews_dates[[1]][ 4 , 1]
as.Date(ews_dates_array[ 3 , 4 , 5 ],origin="1970-01-01")
ews_dates_array[ 3 , 4 , 5 ] == ews_dates[[3]][ 2 , 5]

#' Note one hospitalisation shifted (n-10) date is not transferred to array because there are two dates in the same wave band
#' and one date does not have a time lead/lag because it is in the wave 1 band and there is no start date for wave 1

#' 2 - Cluster logistic growth rates as leading indicator
#'#' Dim = 32 leading indicators, 9 waves and 7 leading indicator statistics. 
#' Data in the array is the EWS date 
as.Date( ews_dates_array[ 1 , 4 , 1 ] , origin = "1970-01-01" )
ews_dates_array[ 1 , 4 , 1 ] == ews_dates[[ 1 ]][ 2 , 1]
as.Date( ews_dates_array[ 3 , 4 , 5 ] , origin = "1970-01-01" )
ews_dates_array[ 3 , 4 , 5 ] == ews_dates[[ 3 ]][ 2 , 5]
as.Date( ews_dates_array[ 20 , 4 , 3 ] , origin = "1970-01-01" )
ews_dates_array[ 3 , 4 , 3 ] == ews_dates[[ 3 ]][ 2 , 3 ]

#' 3 - Ct values as leading indicator
#'#' Dim = 36 leading indicators, 9 waves and 7 leading indicator statistics. 
#' Data in the array is the EWS date 
as.Date(ews_dates_array[ 1 , 4 , 1 ],origin="1970-01-01")
ews_dates_array[ 1 , 4 , 1 ] == ews_dates[[1]][ 2 , 1]
as.Date(ews_dates_array[ 5 , 7 , 7 ],origin="1970-01-01")
ews_dates_array[ 5 , 7 , 7 ] == ews_dates[[5]][ 1 , 7]
as.Date(ews_dates_array[ 27 , 5 , 3 ],origin="1970-01-01")
ews_dates_array[ 27 , 5 , 3 ] == ews_dates[[27]][ 2 , 3]
#' Note 9 EWS dates not transferred to array (possibly because there are dates 
#' in the same wave band)and an additional EWS date does not have a time 
#' lead/lag (possibly because it is in the wave 1 band and there is no start 
#' date for wave 1)

#' 5 - Behavioural - CoMix Survey as leading indicator
#'#' Dim = 14 leading indicators, 9 waves and 7 leading indicator statistics. 
#' Data in the array is the EWS date 
as.Date(ews_dates_array[ 1 , 2 , 1 ],origin="1970-01-01")
ews_dates_array[ 1 , 2 , 1 ] == ews_dates[[1]][ 1 , 1]
as.Date(ews_dates_array[ 5 , 6 , 3 ],origin="1970-01-01")
ews_dates_array[ 5 , 6 , 3 ] == ews_dates[[5]][ 1 , 3]
as.Date(ews_dates_array[ 4 , 4 , 3 ],origin="1970-01-01")
ews_dates_array[ 4 , 4 , 3 ] == ews_dates[[4]][ 1 , 3]
#' All EWS dates transferred to array and all used to calculate lead/lag times

#' 6 - Google mobility as leading indicator
#'#' Dim = 7 leading indicators, 9 waves and 7 leading indicator statistics. 
#' Data in the array is the EWS date 
as.Date(ews_dates_array[ 1 , 4 , 1 ],origin="1970-01-01")
ews_dates_array[ 1 , 4 , 1 ] == ews_dates[[1]][ 4 , 1]
as.Date(ews_dates_array[ 3 , 2 , 3 ],origin="1970-01-01")
ews_dates_array[ 3 , 2 , 3 ] == ews_dates[[3]][ 1 , 3]
as.Date(ews_dates_array[ 3 , 3 , 3 ],origin="1970-01-01")
ews_dates_array[ 3 , 3 , 3 ] == ews_dates[[3]][ 2 , 3]
#' Note, six EWS dates not transferred to array (possibly because there are
#' two dates in the same wave band) and a further 15 dates do not have a time
#' lead/lag (possibly because they are in the wave 1 band and there is no start
#' date for wave 1)


#' Array Test below applies to all leading indicators
#'Method assumes only one date per wave per statistic per leading indicator, so
#'check that all dates captured in dates array...
total_ews_dates <- 0
for (counter in 1 : length( ews_dates ) ) {
  temp <- sum( !is.na( ews_dates[[ counter ]] ))
  total_ews_dates <- total_ews_dates + temp
}
message( paste( "A total of " , total_ews_dates , " EWS dates have been identified.") )

total_ews_dates_in_array <- sum( !is.na( ews_dates_array[  , , ] ))
message( paste( "A total of " , total_ews_dates_in_array , " EWS dates have been transferred to the array.") )

if ( total_ews_dates == total_ews_dates_in_array ){
  message("All EWS dates have been transferred to the array.")
} else {
  message( paste( "NOT all EWS dates have been transferred to the array. There are "
                  , total_ews_dates - total_ews_dates_in_array, " out of ", total_ews_dates
                  , " dates missing. Check list of leading indicators." ) )
  #for (counter in 1 : length( ews_dates ) ) {
  #  temp <- sum( !is.na( ews_dates[[ counter ]] ))
  #  total_ews_dates <- total_ews_dates + temp
  #}
  
}
#' ...and in diff array
total_ews_in_diff_array <- sum( !is.na( ews_diff_array[  , , ] ))
message( paste( "A total of " , total_ews_in_diff_array , " lead/lag calculations have been made in the 'diff' array.") )

if ( total_ews_dates == total_ews_in_diff_array ){
  message("All EWS dates have been used to calculate a lead/lag relative to the wave start date.")
} else {
  message( paste( "NOT all EWS dates have been used to calculate a lead/lag time. There are "
                  , total_ews_dates - total_ews_in_diff_array, " out of ", total_ews_dates
                  , " dates missing. Check list of leading indicators." ) )
  #for (counter in 1 : length( ews_dates ) ) {
  #  temp <- sum( !is.na( ews_dates[[ counter ]] ))
  #  total_ews_dates <- total_ews_dates + temp
  #}
  
}

################################################################################
#' 2 - Plot data with time on the HORIZONTAL

#' 2.1 - Test - hospitalisation data shifted
lead_ind_start = 1 #' Change to plot desired leading indicators
lead_ind_end = 4 #' Change to plot desired leading indicators

#' 2.2 - Cluster logistic growth rates
lead_ind_start = 1 #1 #8 #17 #24 #' Change to plot desired leading indicators
lead_ind_end = 6 #7 #16 #23 #32 #' Change to plot desired leading indicators

#' 2.3 - Ct values
lead_ind_start = 1 #1 #10 #19 #28 #' Change to plot desired leading indicators
lead_ind_end = 9 #9 #18 #27 #36 #' Change to plot desired leading indicators

#' 2.4 - PCR positivity rates in England (pillar 2) 
lead_ind_start = 1 #' Change to plot desired leading indicators
lead_ind_end = 2 #' Change to plot desired leading indicators

#' 2.5 - Behavioural - CoMix Survey
lead_ind_start = 1 #1 #5  #' Change to plot desired leading indicators
lead_ind_end = 4 #4 #13  #' Change to plot desired leading indicators

#' 2.6 - Behavioural - Google mobility
lead_ind_start = 1 #' Change to plot desired leading indicators
lead_ind_end = 7 #' Change to plot desired leading indicators

#' Plot applies to all leading indicators
waves_to_plot = c(2,3,4,5,6,7,8) #' Change to plot desired waves

par( mfrow = c( 9 , length( waves_to_plot ) ) ) #' 9 leading indicator types on the vertical and 7 waves on the horizontal
par( mar = c( 2 , 2 , 1.5 , 1 ) ) #margins: bottom, left, top, right
plot_colours = c("red","blue","green","brown1","cadetblue","darkviolet","orange") 
#x_limits = c( min( ews_diff_array , na.rm = TRUE ) , max( ews_diff_array , na.rm = TRUE) )
#x_limits = c( -plyr::round_any( max( ews_diff_array , na.rm = TRUE) , 10 , f = ceiling ) , -plyr::round_any( min( ews_diff_array , na.rm = TRUE ) , 10 , f = floor ) )#' Use negative as want to show -ve as leading and +ve as lagging
x_limits = c( -100 , 100 )

for ( lead_ind_i in lead_ind_start : lead_ind_end ){
  for (wave_j in waves_to_plot ){  
    x = -ews_diff_array[ lead_ind_i , wave_j ,  ] #' Use negative as want to show -ve as leading and +ve as lagging
    y = seq(1,7,1)
    #message(lead_ind_i, wave_j)
    
    if ( ( lead_ind_i == lead_ind_start ) && ( wave_j == waves_to_plot[ 1 ] ) ){ #'top left corner plot
      par(xaxt="n", yaxt="s")
      plot( x, y, pch=21, col =  plot_colours, bg = plot_colours, cex=2, xlim = x_limits, main = paste("Wave", wave_j))
      #lines( replicate(7,0) , seq(1,7,1) )
      abline( v = 0 )
      #message("yes")
      
    } else if ( ( lead_ind_i == lead_ind_start ) && ( wave_j != waves_to_plot[ 1 ] ) ){ #'remainder of top row plots
      par(xaxt="n", yaxt="n")
      plot( x, y, pch=21, col =  plot_colours, bg = plot_colours, cex=2, xlim = x_limits, main = paste("Wave", wave_j))
      #lines( seq(1,7,1) , replicate(7,0))
      abline( v = 0 )
      #message("yes")
      
    } else if ( ( lead_ind_i == lead_ind_end ) && (wave_j == waves_to_plot[ 1 ] ) ){ #'bottom left corner plot
      par(xaxt="s", yaxt="s")
      plot( x, y, pch=21, col =  plot_colours, bg = plot_colours, cex=2, xlim = x_limits  )
      #lines( seq(1,7,1) , replicate(7,0))
      abline( v = 0 )
      #message("yes")
      
    } else if ( ( lead_ind_i %in% seq( lead_ind_start + 1 , lead_ind_end - 1 , 1 ) ) && (wave_j == waves_to_plot[ 1 ]) ){ #'left column plots but not top or bottom
      par(xaxt="n", yaxt="s")
      plot( x, y, pch=21, col =  plot_colours, bg = plot_colours, cex=2, xlim = x_limits  )
      #lines( seq(1,7,1) , replicate(7,0))
      abline( v = 0 )
      #message("yes")
      
    } else if ( ( lead_ind_i == lead_ind_end ) && ( wave_j != waves_to_plot[ 1 ] ) ){ #'remainder of bottom row plots
      par(xaxt="s", yaxt="n")
      plot( x, y, pch=21, col = plot_colours, bg = plot_colours, cex=2, xlim = x_limits)
      #lines( seq(1,7,1) , replicate(7,0))
      abline( v = 0 )
      #message("yes")
      
    } else if ( ( lead_ind_i %in% seq( lead_ind_start + 1 , lead_ind_end - 1 , 1 ) ) && (wave_j %in% waves_to_plot[2:length(waves_to_plot)] ) ){ #' inner plots
      par(xaxt="n",yaxt="n")
      plot( x, y, pch=21, col =  plot_colours, bg = plot_colours, cex=2, xlim = x_limits )
      #lines( seq(1,7,1) , replicate(7,0))
      abline( v = 0 )
      #message("yes")
    }
  }
}

################################################################################
#' Plot hospitalisation data with wave start dates and bands for comparison with 
#' EWS dates

#' Load case/hospitalisation data
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs'
#filename <- "data_2022-May-09 - UK Covid-19 cases by sample date.csv"
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
#dat_type <- "cases"
dat_type <- "hospitalisations"
# data_load() @ C:\Users\kdrake\GitHub\Early_Warning_Signal
dat_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)

par(xaxt="s", yaxt="s")
par( mar = c( 5 , 5 , 5 , 5 ) ) #margins: bottom, left, top, right

#' plot hospitalisation data
plot( dat_df$time
      , dat_df$cases
      , typ="l"
      , xlab = "Date"
      , ylab = "UK Covid-19 hospitalisations"
      , col="blue")
#' Add wave bands
for ( r in 1 : nrow( wave_bands ) ){
  abline(v = lubridate::decimal_date( wave_bands$band_start[ r ] ) , lty = 2)
  abline(v = lubridate::decimal_date( wave_bands$band_end[ r ] ) , lty = 2)
}
#' Add wave start dates
sd = lubridate::decimal_date( wave_bands$wave_start )
points(   lubridate::decimal_date( wave_bands$wave_start )
          , dat_df$cases[ match( sd , dat_df$time ) ]
          , col="red"
          , cex = 3
)

