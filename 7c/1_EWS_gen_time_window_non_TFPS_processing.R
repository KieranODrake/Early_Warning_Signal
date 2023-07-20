#' Purpose: Prepare non-TFPS leading indicator data ready for generation of early warning signals (EWS) 
#' (non-TFPS, for TFPS EWS see 'EWS_calc_threshold_lgr_th_XXX_pval_filter_XXX_HPC_array.R')
#' Author: Kieran Drake
#' Date: March 2023

#' 1 - Load data for each leading indicator type and format in dataframes with date and different variations in columns
#' 2 - Compute 'robust' z-scores for all leading indicator data
#' 3 - Compute whether z-score above or below EWS threshold
#' 4 - Load waves definitions: start/inflection date, Rt critical transition dates, wave band dates
#' 5 - For each leading indicator, wave, and EWS threshold, want to know number of TP, FP, TN and FN as well as earliest TP.
library( magrittr )
library( lubridate )
#library( ggplot2 )
library( stringr )
library( zoo )
library( data.table )
library( sqldf )
library( zoo )

###' 1 - Load data
#' 
#' Functions to load
#' data_load() @ C:\Users\kdrake\OneDrive - Imperial College London\Documents\GitHub\Early_Warning_Signal

##' 1(a) Load case/hospitalisation data
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs'
#filename <- "data_2022-May-09 - UK Covid-19 cases by sample date.csv"
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
#dat_type <- "cases"
dat_type <- "hospitalisations"
#' data_load() @ C:\Users\kdrake\GitHub\Early_Warning_Signal
dat_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)

#' Define 'shift' function
shift <- function( x , n ){
  c( x[ -( seq( n ) ) ], rep( NA , n ) )
}
hosp = data.frame( "date" = dat_df$date 
                 , "hosp_shift_00" = dat_df$cases
                 , "hosp_shift_10" = shift(dat_df$cases, 10)
                 , "hosp_shift_20" = shift(dat_df$cases, 20) 
                 , "hosp_shift_30" = shift(dat_df$cases, 30) 
)

##' 1(b) load SARS-CoV-2 Ct values (PCR cycle threshold)
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files')
Ct_p2_mean_df   = fread( "Ct_p2_mean_df_v2.csv" )   ;  Ct_p2_mean_df[,1]   = NULL ; Ct_p2_mean_df   = data.frame( Ct_p2_mean_df )
Ct_p2_median_df = fread( "Ct_p2_median_df_v2.csv" ) ;  Ct_p2_median_df[,1] = NULL ; Ct_p2_median_df = data.frame( Ct_p2_median_df )

##' 1(c) load SARS-CoV-2 positivity rates
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/Positivity rates/Outputs/")
positivity_rate_df = readRDS( "sc2_positivity_rate_England.rds" ) ; positivity_rate_df = positivity_rate_df[,c(1,4,5)]
  
##' 1(d) load Behavioural - CoMix survey
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/UK CoMix data 2022/CoMix data')
ews_base = fread( "2022-03-02_bs_means_2w_open.csv" )
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
comix_df = ews_base
rm(ews_00_04,ews_05_11,ews_05_17,ews_18_29,ews_18_59,ews_30_39,ews_40_49,ews_50_59,ews_60_69,ews_60_plus,ews_70_plus,ews_All,ews_All_adults,ews_base)
#' Because the CoMix data is weekly no leading indicator reset dates are identified
#' Therefore need to interpolate to generate dataset of daily datapoints
#' Create list of dates in between weekly measurements
comix_min_date = min(comix_df[,1]) ; comix_max_date = max(comix_df[,1])
comix_cal_dates = seq(comix_min_date,comix_max_date,1)
'%ni%' = Negate('%in%')
comix_cal_dates_missing = comix_cal_dates[ (comix_cal_dates %ni% comix_df[,1]) ]
comix_cal_dates_missing_df = data.frame(  "date" = comix_cal_dates_missing )
#' Add columns for CoMix measurements
for (i in 2 : ncol( comix_df ) ){
  comix_cal_dates_missing_df[ , i ] <- NA
  names( comix_cal_dates_missing_df )[ i ] <- names( comix_df )[ i ]
}
#' Merge rows with weekly data and empty rows for dates in between
comix_interpolate_df = merge( comix_df , comix_cal_dates_missing_df 
                              #, by = intersect( names( comix_df ) , names( comix_cal_dates_missing_df ) )
                              , by.x = names( comix_df ) , by.y = names( comix_cal_dates_missing_df )
                              , all.x = TRUE, all.y = TRUE)
#' Interpolate values in between weekly data
comix_na_approx = comix_interpolate_df
#comix_spline    = comix_interpolate_df
for (i  in 2 : ncol( comix_interpolate_df ) ){
  message(i)
  row_first = which( !is.na( comix_interpolate_df[,i] ) )[ 1 ] #' Find which row is the first to have a value
  row_last = max( which( !is.na( comix_interpolate_df[,i] ) ) )
  comix_na_approx[ row_first:row_last , i ] = na.approx( comix_interpolate_df[ row_first:row_last , i ] )
#  comix_spline[    row_first:row_last , i ] = spline( comix_interpolate_df[row_first:row_last,1] ,   comix_interpolate_df[ row_first:row_last , i ] )
}
#' Check interpolation
plot( comix_na_approx[,1],comix_na_approx[,2])
points( comix_df[,1],comix_df[,2],col="red",pch=16)
#' Remove intermediate dataframes
rm(comix_cal_dates_missing_df,comix_interpolate_df,comix_max_date,comix_min_date,comix_cal_dates,comix_cal_dates_missing)

##' 1(e) load Behavioural - Google mobility
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/Google Mobility')
ews_base = fread( "Google mobility data - UK to 12 June 2022.csv" )
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
google_mobility_df <- ews_base ; rm( ews_base )
#plot(google_mobility_df[,1],google_mobility_df[,2])
#' Add 7-day rolling mean for Google mobility data
for(i in 2:8){
  google_mobility_df[ , i+7 ] =  data.table::frollmean( google_mobility_df[ , i ] , n = 7 , fill = NA )
  names( google_mobility_df )[ i+7 ] <- paste0( names( google_mobility_df )[ i ],"_7d_mean")
}

#' Ensure date column name is the same for all dataframes as will be used for merge
names(hosp)[1] <- "date" ; names(positivity_rate_df)[1] <- "date" ; names(Ct_p2_mean_df)[1] <- "date" ; names(Ct_p2_median_df)[1] <- "date" ; names(comix_na_approx)[1] <- "date" ; names(google_mobility_df)[1] <- "date" 

#' Shift data to account for lag in publication of leading indicator data
lag_hosp       = 7 #' Estimate
lag_positivity = 7
lag_ct         = 7
lag_comix      = 7+3 #' The time series uses the mid-date of a 2-week long survey period and the results were typically published 3 days after the end date of the survey. Report 5 https://cmmid.github.io/topics/covid19/reports/comix/CoMix%20Weekly%20Report%205.pdf 
lag_google     = 3 #' Takes 2-3 days to produce the data https://www.google.com/covid19/mobility/data_documentation.html?hl=en

hosp$date = hosp$date + lag_hosp
positivity_rate_df$date = positivity_rate_df$date + lag_positivity
Ct_p2_mean_df$date = Ct_p2_mean_df$date + lag_ct
Ct_p2_median_df$date = Ct_p2_median_df$date + lag_ct
comix_na_approx$date = comix_na_approx$date + lag_comix
google_mobility_df$date = google_mobility_df$date + lag_google

##' Trim all non-genomic leading indicators to the maximum time period which all leading indicators cover
date_min = max( as.Date("2020-08-14") , min(comix_na_approx$date) , min(Ct_p2_mean_df$date) , min(Ct_p2_median_df$date) , min(google_mobility_df$date), min(hosp$date) , min(positivity_rate_df$date) )
date_max = min( as.Date("2022-03-29") , max(comix_na_approx$date) , max(Ct_p2_mean_df$date) , max(Ct_p2_median_df$date) , max(google_mobility_df$date), max(hosp$date) , max(positivity_rate_df$date) )

trim_comix_na_approx = subset( comix_na_approx    , ( ( comix_na_approx$date    >= date_min ) & ( comix_na_approx$date    <= date_max ) ) )
trim_ct_p2_mean      = subset( Ct_p2_mean_df      , ( ( Ct_p2_mean_df$date      >= date_min ) & ( Ct_p2_mean_df$date      <= date_max ) ) )
trim_ct_p2_median    = subset( Ct_p2_median_df    , ( ( Ct_p2_median_df$date    >= date_min ) & ( Ct_p2_median_df$date    <= date_max ) ) )
trim_hosp            = subset( hosp               , ( ( hosp$date               >= date_min ) & ( hosp$date               <= date_max ) ) )
trim_google_mobility = subset( google_mobility_df , ( ( google_mobility_df$date >= date_min ) & ( google_mobility_df$date <= date_max ) ) )
trim_positivity_rate = subset( positivity_rate_df , ( ( positivity_rate_df$date >= date_min ) & ( positivity_rate_df$date <= date_max ) ) )

#' Merge leading indicator data sets on 'date' column
li_list <- list( trim_hosp , trim_positivity_rate , trim_ct_p2_mean , trim_ct_p2_median , trim_comix_na_approx , trim_google_mobility )
li_df <- li_list[[1]]
for (i in 2:length( li_list ) ) {
  li_df <- merge( li_df ,  li_list[[i]], by="date")
}

#' Replace Inf values with NA and interpolate to replace all NAs 
for(i in 2 : ncol( li_df ) ){
  inf_index = which( is.infinite( li_df[ , i ] ) )
  li_df[ inf_index , i ] <- NA
  li_df[ , i ] = na.approx( li_df[ , i ] )
}


##' 2 - Compute 'robust' z score for all non-genomic leading indicators
#' Loop through time series for each set of variables in the leading indicator file
li_z_df = li_df
li_z_df[ , 2 : ncol( li_df ) ] = NA
for ( var_type in 2 : ncol( li_df ) ){
  #' Process EWS time series (except test case - shifted hospitalisations)
  ews = li_df[ , c( 1 , var_type ) ] # Select date column (1) and variable column (var_type)
  colnames( ews ) <- c( "date" , "cases" ) # Rename columns
  #ews = ews[ !is.na( ews$cases ) , ]# Remove days with NA values
  ews$date <- as.Date( ews$date , format = "%d/%m/%Y") # Change format of date
  ews$time <- lubridate::decimal_date( ews$date ) # Add column for decimal date
  ews$wday <- lubridate::wday( ews$date ) # Add column for day of week
  #ews = ews[ Reduce( '&' , lapply( ews , is.finite ) ) , ] # Remove days with Inf values
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
  li_z_df[ , var_type ] <- ews$running_z_score_robust
  message( names(li_df )[ var_type ] )
}   
plot(li_z_df[,1],rowMeans(li_z_df[,2:65]),ylim=c(-5,+5),typ="l")
abline(h=0)
par( mfrow=c(1,1))
plot(li_z_df[,1],li_z_df[,2],ylim=c(-50,+30),typ="l")
for(i in 3:72){
  lines(li_z_df[,1],li_z_df[,i],col=rainbow(72)[i])
  #invisible( readline( prompt="Press [enter] to continue" ) )
}
#' Closer look at BA.2
plot(li_z_df[,1],li_z_df[,2],ylim=c(-10,+10),xlim=c(wave_bands_start[7],wave_bands_end[7]), typ="l")
for(i in 3:72){
  lines(li_z_df[,1],li_z_df[,i],col=rainbow(72)[i])
  #invisible( readline( prompt="Press [enter] to continue" ) )
}
abline(v=wave_start_dates[7] )
abline(v=c(wave_start_dates[7]+true_ews_start , wave_start_dates[7] + true_ews_end ), lty=2)
              
       
for(i in 2:72){
  plot(li_z_df[,1],li_z_df[,i],col=rainbow(72)[i],ylim=c(-50,30),main=names(li_z_df)[i],typ="l",lwd=2)
  lines(li_df[,1],li_df[,i],col=rainbow(72)[i+5],ylim=c(-50,30),main=names(li_z_df)[i],typ="l")
  abline(h=0)
  #invisible( readline( prompt="Press [enter] to continue" ) )
}


plot(li_z_df[,1],li_z_df[,2],ylim=c(-30,+30),typ="l")
abline(h=0)

#' 3 - Compute whether z-score above or below EWS threshold
#' Create an array. 71 columns for different leading indicators, 563 rows for dates and 101 z-dim for EWS threshold levels
#' Perhaps do more threshold levels? -5 to +5? instead of just 0 to +5?
li_z_ews_array = array(data=NA
                       ,dim=c( nrow(li_z_df) #' dates
                               ,ncol(li_z_df) - 1 #' leading indicator 'robust' z-score values -ignore first column which is the dates
                               ,length(seq(-5,+5,0.05)) #' Whether above or below EWS threshold
                       )
)
ews_count = 0
for( ews_th in seq(-5,5,0.05) ){
  ews_count = ews_count + 1 
  temp = !(li_z_df[,2:ncol(li_z_df)] <= ews_th ) #' TRUE if above EWS threshold and FALSE if below
  #' Replace Inf and NA values with FALSE 
  for(i in 1 : ncol( temp ) ){
    inf_na_index = which( is.infinite( temp[ , i ] ) | is.na( temp[ , i ] ) )
    temp[ inf_na_index , i ] <- FALSE
  }
  li_z_ews_array[ , , ews_count ] = temp
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

Rt_critical_transitions <- c(                      #' 1 Wuhan
  as.Date("2020-09-06") #' 2 B.1.177
  ,as.Date("2020-12-13") #' 3 Alpha
  ,as.Date("2021-05-22") #' 4 Delta
  ,as.Date("2021-08-19") #' 5 Delta
  ,as.Date("2021-10-13") #' 6 Delta
  ,as.Date("2021-11-25") #' 7 Omicron
  ,as.Date("2022-03-11") #' 8 Omicron
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
#' Check dates
plot(x=seq(wave_bands_start[1],wave_bands_end[7],1)
     ,y=replicate( length( seq(wave_bands_start[1],wave_bands_end[7],1)),0))
for (i in 1:length(wave_bands_start)){
  abline(v=wave_bands_start[i],lty=2)
  abline(v=wave_bands_end[i],lty=2)
  abline(v=wave_start_dates[i],lty=1,col="red")
  abline(v=wave_start_dates[i]+true_ews_start,lty=2,col="red")
  abline(v=wave_start_dates[i]+true_ews_end,lty=2,col="red")
  abline(v=Rt_critical_transitions[i],lty=1,col="orange")
  abline(v=Rt_critical_transitions[i]+true_ews_start,lty=2,col="orange")
  abline(v=Rt_critical_transitions[i]+true_ews_end,lty=2,col="orange")
  invisible( readline( prompt="Press [enter] to continue" ) )
}
  

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
true_ews_end = + 5 #5 #' t+5 as per Proverbio et al 
true_start_dates = list() ; true_rt_dates = list()
for(i in 1:7){
  #' Dates where value above threshold would be a TRUE Positive (TP) EWS and below threshold would be a FALSE negative (FN) EWS
  true_start_dates[[i]] = seq( wave_start_dates[i] + true_ews_start , wave_start_dates[i] + true_ews_end , 1 )
  message(wave_start_dates[i] + true_ews_start ," ", wave_start_dates[i] + true_ews_end)
  true_rt_dates[[i]] = seq( Rt_critical_transitions[i] + true_ews_start , Rt_critical_transitions[i] + true_ews_end , 1 )
}
#' On dates other than those defined above, the classification is reversed:
#' False positive (FP) if above EWS threshold and TRUE negative (TN) if below EWS threshold
#' Apply to array where TRUE simply indicates above the EWS threshold and FALSE indicates below or equal to EWS threshold
true_start_dates_index = which( li_z_df[,1] %in% unlist(true_start_dates) )
true_rt_dates_index = which( li_z_df[,1] %in% unlist(true_rt_dates) )

li_z_ews_sd_classification_array = li_z_ews_array #' Initialise new array for classifying signals using existing array as template
li_z_ews_rt_classification_array = li_z_ews_array #' Initialise new array for classifying signals using existing array as template

for( k in 1:nrow(li_z_df) ){ #' loop through dates
  for( m in 1:(ncol(li_z_df)-1) ){ #' loop through leading indicators
    for( n in 1:dim(li_z_ews_array)[3] ){ #' loop through EWS thresholds
      #' Classify based on wave start dates
      if( (k %in% true_start_dates_index ) & ( li_z_ews_array[k,m,n] == TRUE ) ){ #' WITHIN EWS wave range & ABOVE EWS threshold
        li_z_ews_sd_classification_array[k,m,n] = "TP"
      }
      if( (k %in% true_start_dates_index ) & ( li_z_ews_array[k,m,n] == FALSE ) ){ #' WITHIN EWS wave range & BELOW EWS threshold
        li_z_ews_sd_classification_array[k,m,n] = "FN"
      }
      if( (k %ni% true_start_dates_index ) & ( li_z_ews_array[k,m,n] == TRUE ) ){ #' OUTSIDE EWS wave range & ABOVE EWS threshold
        li_z_ews_sd_classification_array[k,m,n] = "FP"
      }
      if( (k %ni% true_start_dates_index ) & (li_z_ews_array[k,m,n] == FALSE ) ){ #' OUTSIDE EWS wave range & BELOW EWS threshold
        li_z_ews_sd_classification_array[k,m,n] = "TN"
      }
      #' Classify based on Rt critical transition dates
      if( (k %in% true_rt_dates_index ) & ( li_z_ews_array[k,m,n] == TRUE ) ){ #' WITHIN EWS wave range & ABOVE EWS threshold
        li_z_ews_rt_classification_array[k,m,n] = "TP"
      }
      if( (k %in% true_rt_dates_index ) & ( li_z_ews_array[k,m,n] == FALSE ) ){ #' WITHIN EWS wave range & BELOW EWS threshold
        li_z_ews_rt_classification_array[k,m,n] = "FN"
      }
      if( (k %ni% true_rt_dates_index ) & ( li_z_ews_array[k,m,n] == TRUE ) ){ #' OUTSIDE EWS wave range & ABOVE EWS threshold
        li_z_ews_rt_classification_array[k,m,n] = "FP"
      }
      if( (k %ni% true_rt_dates_index ) & (li_z_ews_array[k,m,n] == FALSE ) ){ #' OUTSIDE EWS wave range & BELOW EWS threshold
        li_z_ews_rt_classification_array[k,m,n] = "TN"
      }
    }
  }
}

#' 5 - For each leading indicator, wave, and EWS threshold, want to know number of TP, FP, TN and FN as well as earliest TP.
#'  classification types
#' 201 EWS thresholds
#' 7 waves + 1 total + 1 total ex.B.1.177 + 1 total ex.B.1.177 & BA.2

#' Create range of dates associated with each wave
date_index_list = list()
for(p in 1:7){
  date_index_temp = which( li_z_df[,1] %in% seq( wave_bands_start[ p ], wave_bands_end[ p ] , 1 ) )
  assign( paste0( "date_index_w" , p+1 ), date_index_temp )
  date_index_list[[p]] <- date_index_temp
  names(date_index_list)[p] <- paste0( "date_index_w" , p+1 )
}
#' Create date range for total
date_index_total            <- which( li_z_df[,1] %in% seq( wave_bands_start[ 1 ], wave_bands_end[ 7 ] , 1 ) )
#' Create date range for total ex. B.1.177 (not enough data)
date_index_total_ex_B.1.177 <- which( li_z_df[,1] %in% seq( wave_bands_start[ 2 ], wave_bands_end[ 7 ] , 1 ) )
#' Create date range for total ex. B.1.177 and BA.2 (not enough data)
date_index_total_ex_B.1.177_BA.1 <- which( li_z_df[,1] %in% seq( wave_bands_start[ 2 ], wave_bands_end[ 6 ] , 1 ) )
#' Add to list of lists
date_index_list[[ 8]] = date_index_total                 ; names(date_index_list)[ 8] <- "date_index_total"
date_index_list[[ 9]] = date_index_total_ex_B.1.177      ; names(date_index_list)[ 9] <- "date_index_total_ex_B.1.177"
date_index_list[[10]] = date_index_total_ex_B.1.177_BA.1 ; names(date_index_list)[10] <- "date_index_total_ex_B.1.177_BA.1"
