#' Purpose: Calculate early warning signals (EWS) using leading indicator data.
#' A range of leading indicator thresholds are used for generating EWS. 
#' Absolute value of the leading indicator data is used.
#' EWS are assessed as True Positive or False Positive. This is  based on whether
#' the dominant variant during the proximal (this definition may need some 
#' tightening - currently dates between c.midway down RHS of previous peak and
#' midway down RHS of peak in question) infection wave matches any of:
#'  ~ variants listed in top 5 most frequent in terms of the number of clusters where it is the most frequent
#'  ~  most frequent in clusters with 5 largest logistic growth rates. 
#' Author: Kieran Drake
#' Date: October/November 2022
#' 

library( magrittr ) 
library( lubridate ) 
#library( ggplot2 )
library( stringr )
library( zoo )
library( data.table )
library( sqldf )

#' array number from High Performance Computing (HPC) array job number
params_run <- commandArgs( trailingOnly = TRUE )
#' values to be varied in HPC array run
parent_sub_lgr_threshold = 0.60 #' 0.65 0.70 0.75 0.80 0.85 0.90 0.95
p_val_thresholds = 0.01 #10000 #, 0.05 , 0.01 


#########' 1 - Load data ####################
#' Wave start dates (average where a wave has more than one start date as 
#' produced from the 12 models incorporating 'low resolution wave filter' - see
#' 21 June 2022 update slides)
wave_start_dates <- c(            #' 1 Wuhan
                       as.Date("2020-08-19") #' 2 B.1.177
                       ,as.Date("2020-11-29") #' 3 Alpha
                       ,as.Date("2021-05-11") #' 4 Delta
                       ,as.Date("2021-08-03") #' 5 Delta
                       ,as.Date("2021-09-27") #' 6 Delta
                       ,as.Date("2021-11-26") #' 7 Omicron
                       ,as.Date("2022-02-21") #' 8 Omicron
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

#' Initialise Early Warning Signal (EWS) results dataframe
ews_results_df = data.frame(  "TFPS_cluster_min_age" = NA
                              , "TFPS_cluster_max_age" = NA
                              , "TFPS_cluster_min_descendants" = NA
                              , "LGR_p_value_threshold" = NA
                              , "parent_sub_lgr_threshold" = NA
                              , "leading_indicator_type" = NA
                              , "leading_indicator_value" = NA
                              , "EWS_threshold" = NA
                              , "EWS_date" = NA
)
new_row = ews_results_df[ 1 , ]
#' Initialise format for final results dataframe combining the other dataframes
results_df = data.frame()
#results_df = data.frame(  "wave_n" = NA
#                        , "wave_name" = NA
#                        , "pango" = NA
#                        , "wave_inflection_start_date" = NA
#                        , "wave_band_start" = NA 
#                        , "wave_band_end" = NA
#                        , "TFPS_cluster_min_age" = NA
#                        , "TFPS_cluster_max_age" = NA
#                        , "TFPS_cluster_min_descendants" = NA
#                        , "LGR_p_value_threshold" = NA
#                        , "parent_sub_lgr_threshold" = NA
#                        , "leading_indicator_type" = NA
#                        , "leading_indicator_value" = NA
#                        , "EWS_threshold" = NA
#                        , "EWS_date" = NA
#                        , "lead_lag_days" = NA
#                        , "highest_freq_variant_1" = NA, "n_clusters_1" = NA
#                        , "highest_freq_variant_2" = NA, "n_clusters_2" = NA
#                        , "highest_freq_variant_3" = NA, "n_clusters_3" = NA
#                        , "highest_freq_variant_4" = NA, "n_clusters_4" = NA
#                        , "highest_freq_variant_5" = NA, "n_clusters_5" = NA
#                        , "variant_1_LGR"=NA, "variant_1_LGR_freq"=NA,"variant_1_LGR_value"=NA,"variant_1_simple_LGR"=NA,"variant_1_simple_LGR_freq"=NA,"variant_1_simple_LGR_value"=NA,"variant_1_GAM_LGR"=NA,"variant_1_GAM_LGR_freq"=NA,"variant_1_GAM_LGR_value"=NA #' Information on variant with highest relevant growth rate
#                        , "variant_2_LGR"=NA,"variant_2_LGR_freq"=NA,"variant_2_LGR_value"=NA,"variant_2_simple_LGR"=NA,"variant_2_simple_LGR_freq"=NA,"variant_2_simple_LGR_value"=NA,"variant_2_GAM_LGR"=NA,"variant_2_GAM_LGR_freq"=NA,"variant_2_GAM_LGR_value"=NA
#                        , "variant_3_LGR"=NA,"variant_3_LGR_freq"=NA,"variant_3_LGR_value"=NA,"variant_3_simple_LGR"=NA,"variant_3_simple_LGR_freq"=NA,"variant_3_simple_LGR_value"=NA,"variant_3_GAM_LGR"=NA,"variant_3_GAM_LGR_freq"=NA,"variant_3_GAM_LGR_value"=NA
#                        , "variant_4_LGR"=NA,"variant_4_LGR_freq"=NA,"variant_4_LGR_value"=NA,"variant_4_simple_LGR"=NA,"variant_4_simple_LGR_freq"=NA,"variant_4_simple_LGR_value"=NA,"variant_4_GAM_LGR"=NA,"variant_4_GAM_LGR_freq"=NA,"variant_4_GAM_LGR_value"=NA
#                        , "variant_5_LGR"=NA,"variant_5_LGR_freq"=NA,"variant_5_LGR_value"=NA,"variant_5_simple_LGR"=NA,"variant_5_simple_LGR_freq"=NA,"variant_5_simple_LGR_value"=NA,"variant_5_GAM_LGR"=NA,"variant_5_GAM_LGR_freq"=NA,"variant_5_GAM_LGR_value"=NA
#                        , "positive_EWS" #' True of false positive depending on whether the variant causing the wave appears in top variants for five clusters with largest growth values 
#                        )

############################################
#' Read in dataframe containing information on highest growth variants and most frequent variants
#' for each set of scan variables (min age, max age, min descendants) and filter variable (LGR p-value) 
#setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_11/analysis")
top_variants_df = readRDS( "top_variants_df.rds" ) ; top_variants_df$tree_date = as.Date(top_variants_df$tree_date ) #readRDS( "top_variants_positive_large_df.rds" )

#colnames( top_variants_df[ 6 ] ) <- c("parent_sub_lgr_threshold")
#colnames ( top_variants_df ) [ colnames ( top_variants_df ) == "int_cluster_LGR_threshold" ] <- "parent_sub_lgr_threshold"

###################
#' Read in leading indicator data
#folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/analysis/p_val_filter_005/dataframes - statistics'
filenames <- c(   "tfps_vlgr_samp.csv"      , "tfps_vlgr_simple_samp.csv"    , "tfps_vlgr_gam_samp.csv" 
                  , "tfps_vlgr_pop.csv"     , "tfps_vlgr_simple_pop.csv"     , "tfps_vlgr_gam_pop.csv"
                  , "tfps_vlgr_wtd.csv"
                  , "tfps_lgr_max.csv"      , "tfps_lgr_simple_max.csv"      , "tfps_lgr_gam_max.csv"
                  , "tfps_lgr_mean.csv"     , "tfps_lgr_simple_mean.csv"     , "tfps_lgr_gam_mean.csv"
                  , "tfps_lgr_wtd_mean.csv" , "tfps_lgr_simple_wtd_mean.csv" , "tfps_lgr_gam_wtd_mean.csv"
                  , "tfps_clock_outlier_max.csv" , "tfps_clock_outlier_mean.csv"
            		  , "tfps_pango_lineage_max_lgr.csv"
)
filenames_prefix <- paste( data.frame( str_split( filenames , pattern = ".csv" ) )[ 1 , ] , "_" , sep = "" )

####################
#'Set parameters for Early Warning Signal (EWS) calculation
#' list of thresholds 
ews_thresholds = seq( 0.00 , 5.00 , 0.05 ) #' this range of EWS thresholds might not be suitable for all leading indicators 
#' Number of consecutive data points required to be above threshold before EWS recorded
n_data_points = 2

#' Logistic growth rate p-value thresholds
#p_val_thresholds = c( 10000 , 0.05 , 0.01 ) 
#' Logistic growth rate (LGR) threshold for replacing external sub-cluster with internal parent cluster
#'**REMEMBER TO ALSO LOOK AT THE SET OF LEADING INDICATORS WHERE NO REPLACEMENTS WERE ATTEMPTED i.e. 'large_cluster_adjust' = FALSE in
#parent_sub_lgr_threshold = 0.60 #' This should be varied

####################
#' Fill ews_results_df dataframe
#start_time = Sys.time()
#' Loop through LGR p-value thresholds
#for (h in 1 : length( p_val_thresholds ) ){
#  if( p_val_thresholds[ h ] == 10000){ folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_10/analysis/p_val_filter_no/dataframes - statistics'}
#  if( p_val_thresholds[ h ] == 0.05 ){ folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_10/analysis/p_val_filter_005/dataframes - statistics'}
#  if( p_val_thresholds[ h ] == 0.01 ){ folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_10/analysis/p_val_filter_001/dataframes - statistics'}

  #' Loop through leading indicator files
#  for ( i in 1 : length( filenames ) ){
    i = params_run
message(filenames)
message(filenames_prefix)
message(i)
    filename = filenames[ as.numeric( i ) ]
    filename_prefix = filenames_prefix[ as.numeric( i ) ]
message(filename," ",filenames[as.numeric(i)])
message(filename_prefix," ",filenames_prefix[ as.numeric(i)])

    #setwd( folder )
    ews_base = data.table::fread( filename ) ;  ews_base[ , 1 ] = NULL ; ews_base = data.frame( ews_base )
    
    #' Loop through time series for each set of variables in the leading indicator file
    for ( var_type in 2 : ncol( ews_base ) ){
      #' Process EWS time series (except test case - shifted hospitalisations)
      ews = ews_base[ , c( 1 , var_type ) ] # Select date column (1) and variable column (var_type)
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
      colnames( ews[ 5 ] ) <- c( "running_z_score" )
      colnames( ews[ 6 ] ) <- c( "running_z_score_robust" )
      for( j in 1 : length( ews_thresholds ) ){
        ews_signal = data.frame()
        #plot( ews[,1] , ews[,6], ylim=c(0,5),typ="l")
        #abline(h=0) ; abline(h=ews_thresholds[j])
        for ( k in n_data_points : dim( ews )[ 1 ] ){
          #' Data set containing the previous n data points required to be above the threshold in order for an EWS to be recorded
          ews_test_dat = ews[ (k - n_data_points + 1) : k , 6 ]
          #' Set first n ews_signal values to 0, as not possible to generate a signal until have at least n data points
          ews_signal[ 1:k , 1 ] = 0
#message("1 ", sum( as.numeric( !is.na( ews_test_dat ) ) ) == n_data_points )
          if ( sum( as.numeric( !is.na( ews_test_dat ) ) ) == n_data_points ){ #' Test if all data points in range are not NA
#message("2 ", sum( as.numeric( ews_test_dat > ews_thresholds[ j ] ) ) == n_data_points )
            if( sum( as.numeric( ews_test_dat > ews_thresholds[ j ] ) ) == n_data_points ){ #' Test if all data points are greater than the EWS threshold
              #'Record an EWS in ews_results_df
              ews_signal[ k , ] = 1
              
              #points(ews[k,1],ews[k,6],col="red",pch=16)
              
              new_row = ews_results_df[ 1 , ]
              new_row$TFPS_cluster_min_age = substr( strsplit( colnames( ews_base[var_type] ) , "_" )[[1]][1] , 5 , 6 )
              new_row$TFPS_cluster_max_age = substr( strsplit( colnames( ews_base[var_type] ) , "_" )[[1]][2] , 5 , 6 )
                md = substr( strsplit( colnames( ews_base[var_type] ) , "_" )[[1]][3] , 3 , nchar( strsplit( colnames( ews_base[var_type] ) , "_" )[[1]][3]) )
#message("3 ", md == "perc" )
                md = ifelse( md == "perc" , "proportional" , md )
              new_row$TFPS_cluster_min_descendants = md
              new_row$LGR_p_value_threshold = p_val_thresholds #[ h ]
              new_row$parent_sub_lgr_threshold = parent_sub_lgr_threshold
              new_row$leading_indicator_type = substr( filename_prefix , 1 , nchar(filename_prefix)-1 )
              new_row$leading_indicator_value = ews[ k , 6 ] #' This is the robust z-score value so should be comparable across other leading indicators
              new_row$EWS_threshold = ews_thresholds[ j ]
              new_row$EWS_date = as.Date( ews$date[ k ] , origin = "1970-01-01" )
              
              ews_results_df = rbind( ews_results_df , new_row )
              
            } else { ews_signal[ k , ] = 0 }
          }
        }
        #message( filename_prefix , " " , p_val_thresholds[ h ] ," ",new_row$TFPS_cluster_min_age," ",new_row$TFPS_cluster_max_age," ",new_row$TFPS_cluster_min_descendants," ", ews_thresholds[ j ] )
      }
    }  
#  }
#}
#end_time = Sys.time()
#message( total_time = end_time - start_time )

#' Remove 1st row of NAs
ews_results_df = ews_results_df[ -1 , ]
#' Transform 'TFPS_cluster_min_descendants' to numeric so matches format in top_variants_df
ews_results_df$TFPS_cluster_min_descendants = as.numeric(ews_results_df$TFPS_cluster_min_descendants)
#' Line above converts TFPS_cluster_min_descendants = 'proportional' to NA so these need to be changed back
for ( m in 1: nrow(ews_results_df) ){
#message("4 ", is.na( ews_results_df$TFPS_cluster_min_descendants[ m ] ) )
  if( is.na( ews_results_df$TFPS_cluster_min_descendants[ m ] ) ){
    ews_results_df$TFPS_cluster_min_descendants[ m ] = "proportional"
  }
}

ews_results_df$EWS_date = as.Date( ews_results_df$EWS_date, origin="1970-01-01" )

#' save dataframe
#setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_10/analysis")
saveRDS( ews_results_df , file = paste0( "ews_results_df_lgr_th_",parent_sub_lgr_threshold,"_pval_th_",p_val_thresholds,"_li_",filenames_prefix[ as.numeric( i ) ],".rds" ) )
#**WILL BE TOO BIG TO SAVE A .CSV**
#write.csv( ews_results_df , file = "ews_results_df.csv") 

#############

#' question about whether ignore further EWS for same wave - probably don't want to ignore them at this stage
#' We would know at the time if the following EWS are for the same wave because they will contain the same variant.

#############
#' Need to sort the EWS by wave and merge with wave_df and top_variants_df

#' Merge ews_results_df and top_variants_df by EWS_date and tree_date respectively but also taking account of the scan variables
library(sqldf)
#results_df = sqldf( "SELECT * 
#                          FROM ews_results_df
#                          LEFT JOIN top_variants_df 
#                          ON ews_results_df.EWS_date = top_variants_df.tree_date 
#                            AND ews_results_df.TFPS_cluster_min_age = top_variants_df.TFPS_cluster_min_age 
#                            AND ews_results_df.TFPS_cluster_max_age = top_variants_df.TFPS_cluster_max_age 
#                            AND ews_results_df.TFPS_cluster_min_descendants = top_variants_df.TFPS_cluster_min_descendants
#                            AND ews_results_df.LGR_p_value_threshold = top_variants_df.LGR_p_value_threshold
#                            AND ews_results_df.parent_sub_lgr_threshold = top_variants_df.parent_sub_lgr_threshold ; " )

names(ews_results_df)[5] <- "parent_sub_LGR_threshold"
#typeof(ews_results_df$parent_sub_LGR_threshold)
ews_results_df$parent_sub_LGR_threshold = as.double(ews_results_df$parent_sub_LGR_threshold)
#typeof(ews_results_df$parent_sub_LGR_threshold)
#typeof(top_variants_df$parent_sub_LGR_threshold)

#results_df = sqldf( "SELECT
#                        ews_results_df.TFPS_cluster_min_age ,
#                        ews_results_df.TFPS_cluster_max_age ,
#                        ews_results_df.TFPS_cluster_min_descendants ,
#                        ews_results_df.LGR_p_value_threshold ,
#                        ews_results_df.parent_sub_LGR_threshold ,
#                        ews_results_df.leading_indicator_type ,
#                        ews_results_df.leading_indicator_value ,
#                        ews_results_df.EWS_threshold ,
#                        ews_results_df.EWS_date ,
#                        top_variants_df.tree_date ,
#                        top_variants_df.highest_freq_variant_1 ,
#                        top_variants_df.n_clusters_1 ,
#                        top_variants_df.highest_freq_variant_2 ,
#                        top_variants_df.n_clusters_2 ,
#                        top_variants_df.highest_freq_variant_3 ,
#                        top_variants_df.n_clusters_3 ,
#                        top_variants_df.highest_freq_variant_4 ,
#                        top_variants_df.n_clusters_4 ,
#                        top_variants_df.highest_freq_variant_5 ,
#                        top_variants_df.n_clusters_5 ,
#                        top_variants_df.variant_1_LGR ,
#                        top_variants_df.variant_1_LGR_freq ,
#                        top_variants_df.variant_1_LGR_value ,
#                        top_variants_df.variant_1_simple_LGR ,
#                        top_variants_df.variant_1_simple_LGR_freq ,
#                       top_variants_df.variant_1_simple_LGR_value ,
#                        top_variants_df.variant_1_GAM_LGR ,
#                        top_variants_df.variant_1_GAM_LGR_freq ,
#                        top_variants_df.variant_1_GAM_LGR_value ,
#                        top_variants_df.variant_2_LGR ,
#                        top_variants_df.variant_2_LGR_freq ,
#                        top_variants_df.variant_2_LGR_value ,
#                        top_variants_df.variant_2_simple_LGR ,
#                        top_variants_df.variant_2_simple_LGR_freq ,
#                        top_variants_df.variant_2_simple_LGR_value ,
#                        top_variants_df.variant_2_GAM_LGR ,
#                        top_variants_df.variant_2_GAM_LGR_freq ,
#                        top_variants_df.variant_2_GAM_LGR_value ,
#                        top_variants_df.variant_3_LGR ,
#                        top_variants_df.variant_3_LGR_freq ,
#                        top_variants_df.variant_3_LGR_value ,
#                        top_variants_df.variant_3_simple_LGR ,
#                        top_variants_df.variant_3_simple_LGR_freq ,
#                        top_variants_df.variant_3_simple_LGR_value ,
#                        top_variants_df.variant_3_GAM_LGR ,
#                        top_variants_df.variant_3_GAM_LGR_freq ,
#                        top_variants_df.variant_3_GAM_LGR_value ,
#                        top_variants_df.variant_4_LGR ,
#                        top_variants_df.variant_4_LGR_freq ,
#                        top_variants_df.variant_4_LGR_value ,
#                        top_variants_df.variant_4_simple_LGR ,
#                        top_variants_df.variant_4_simple_LGR_freq ,
#                        top_variants_df.variant_4_simple_LGR_value ,
#                        top_variants_df.variant_4_GAM_LGR ,
#                        top_variants_df.variant_4_GAM_LGR_freq ,
# #                       top_variants_df.variant_4_GAM_LGR_value ,
#                        top_variants_df.variant_5_LGR ,
#                        top_variants_df.variant_5_LGR_freq ,
#                        top_variants_df.variant_5_LGR_value ,
#                        top_variants_df.variant_5_simple_LGR ,
#                        top_variants_df.variant_5_simple_LGR_freq ,
#                        top_variants_df.variant_5_simple_LGR_value ,
#                        top_variants_df.variant_5_GAM_LGR ,
#                        top_variants_df.variant_5_GAM_LGR_freq ,
#                        top_variants_df.variant_5_GAM_LGR_value
#                      FROM ews_results_df LEFT JOIN top_variants_df 
#                      ON ews_results_df.EWS_date = top_variants_df.tree_date 
#                        AND ews_results_df.TFPS_cluster_min_age = top_variants_df.TFPS_cluster_min_age 
#                        AND ews_results_df.TFPS_cluster_max_age = top_variants_df.TFPS_cluster_max_age 
#                        AND ews_results_df.TFPS_cluster_min_descendants = top_variants_df.TFPS_cluster_min_descendants
#                        AND ews_results_df.LGR_p_value_threshold = top_variants_df.LGR_p_value_threshold 
#                        AND ews_results_df.parent_sub_LGR_threshold = top_variants_df.parent_sub_LGR_threshold; " )

#top_variants_df_filtered = subset( top_variants_df , top_variants_df$parent_sub_LGR_threshold == parent_sub_lgr_threshold & top_variants_df$LGR_p_value_threshold==p_val_thresholds)

#' Due to formatting of the parent_sub_lgr_threshold column (I think) the sqldf merger doesn't work and 
#' so the workaround is to subset top_variants_df first and then join
top_variants_df_filtered = subset( top_variants_df , (top_variants_df$parent_sub_LGR_threshold == parent_sub_lgr_threshold) & (top_variants_df$LGR_p_value_threshold==p_val_thresholds))
merged_temp = sqldf("SELECT * 
                    FROM ews_results_df LEFT JOIN top_variants_df_filtered 
                    ON  ews_results_df.EWS_date = top_variants_df_filtered.tree_date
                    AND ews_results_df.TFPS_cluster_min_age = top_variants_df_filtered.TFPS_cluster_min_age
                    AND ews_results_df.TFPS_cluster_max_age = top_variants_df_filtered.TFPS_cluster_max_age
                    AND ews_results_df.TFPS_cluster_min_descendants = top_variants_df_filtered.TFPS_cluster_min_descendants")
merged_temp[,c(11,12,13,14,15)] = NULL
results_df = merged_temp

#' Add wave data
wave_results_df = sqldf( "SELECT * 
                     FROM results_df
                     LEFT JOIN wave_df 
                     ON EWS_date BETWEEN wave_band_start AND wave_band_end")

#' Calculate lead(-ve)/lag(+ve) time (days) 
lead_lag_days = as.Date( wave_results_df$EWS_date , origin = "1970-01-01" ) - wave_results_df$wave_inflection_start_date
wave_results_df = cbind( wave_results_df , lead_lag_days )

#' Determine whether EWS is true positive or false positive
#' ** Some aliases and/or sub-lineages included ** https://www.pango.network/the-pango-nomenclature-system/statement-of-nomenclature-rules/
#'** Also try filtering for variants where cluster has positive growth rate?** 
positive_EWS = data.frame(NA)
for ( w in 1 : nrow( wave_results_df ) ){
  temp_variant_df = wave_results_df[ w , c("highest_freq_variant_1","highest_freq_variant_2","highest_freq_variant_3","highest_freq_variant_4","highest_freq_variant_5"
                                          , "variant_1_LGR","variant_2_LGR","variant_3_LGR","variant_4_LGR","variant_5_LGR"
                                          , "variant_1_simple_LGR","variant_2_simple_LGR","variant_3_simple_LGR","variant_4_simple_LGR","variant_5_simple_LGR"
                                          , "variant_1_GAM_LGR","variant_2_GAM_LGR","variant_3_GAM_LGR","variant_4_GAM_LGR","variant_5_GAM_LGR")]
  #' This if statement works for lineage and sub-lineage **still need to incorporate aliases**
#message("5 ", sum( grepl( wave_results_df$pango[ w ] , temp_variant_df ) ) > 0 )
#message("6 ", ( wave_results_df$wave_n[ w ] %in% c(4,5,6) ) & ( sum( grepl( "AY.", temp_variant_df ) ) > 0 ) )
#message("7 ", wave_results_df$wave_n[ w ] == 3 & sum( grepl( "Q.", temp_variant_df ) ) > 0 )
  if ( sum( grepl( wave_results_df$pango[ w ] , temp_variant_df ) ) > 0 ){ 
    positive_EWS[ w , 1 ] = TRUE
  #' Below if statement only works for exact match  
  #if ( wave_results_df$pango[ w ] %in% wave_results_df[ w , c("highest_freq_variant_1","highest_freq_variant_2","highest_freq_variant_3","highest_freq_variant_4","highest_freq_variant_5"
  #                                                     , "variant_1_LGR","variant_2_LGR","variant_3_LGR","variant_4_LGR","variant_5_LGR"
  #                                                     , "variant_1_simple_LGR","variant_2_simple_LGR","variant_3_simple_LGR","variant_4_simple_LGR","variant_5_simple_LGR"
  #                                                    , "variant_1_GAM_LGR","variant_2_GAM_LGR","variant_3_GAM_LGR","variant_4_GAM_LGR","variant_5_GAM_LGR")]){ 
  #  positive_EWS[ w , 1 ] = TRUE
  } else if ( ( wave_results_df$wave_n[ w ] %in% c(4,5,6) ) & ( sum( grepl( "AY.", temp_variant_df ) ) > 0 ) ){ #' Also match alias for Delta variant
    positive_EWS[ w , 1 ] = TRUE
  } else if ( wave_results_df$wave_n[ w ] == 3 & sum( grepl( "Q.", temp_variant_df ) ) > 0 ){ #' Also match alias for Alpha variant
      positive_EWS[ w , 1 ] = TRUE
  } else { positive_EWS[ w , 1 ] = FALSE }
}
wave_results_df = cbind( wave_results_df , positive_EWS )
colnames( wave_results_df[ 73 ] ) <- "positive_EWS"
colnames ( wave_results_df ) [ colnames ( wave_results_df ) == "NA." ] <- "positive_EWS"

#' save dataframe
#setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_11/analysis")
saveRDS( wave_results_df , file = paste0( "wave_results_df_lgr_th_",parent_sub_lgr_threshold,"_pval_th_",p_val_thresholds,"_li_",filenames_prefix[ as.numeric( i ) ],".rds" ) )
#**WILL BE TOO BIG TO SAVE A .CSV**
#write.csv( wave_results_df , file = "wave_results_df.csv")

#**Some of these may actually be true positives but a different pango id is being used e.g. B.1.617.2 = AY.4 for Delta**
#**Also should we count sub-lineages like B.1.617.2.1 = AY.4.1 as a true positive?**
#**Also, may need to rethink wave bands in case an early warning signal doesn't fall inside the arbitrary band**
#**but not sure how definitively decide which variant is causing the EWS to then match it to a wave**

###################
#' Tabulate, by wave and variables, the earliest true positive and number of false negatives

#' Initialise dataframe
wave_results_analysis_df = data.frame(  "wave_n" = NA
                                      , "wave_name" = NA
                                      , "pango" = NA
                                      , "wave_inflection_start_date" = NA
                                      , "wave_band_start" = NA 
                                      , "wave_band_end" = NA
                                      , "TFPS_cluster_min_age" = NA
                                      , "TFPS_cluster_max_age" = NA
                                      , "TFPS_cluster_min_descendants" = NA
                                      , "LGR_p_value_threshold" = NA
                                      , "parent_sub_LGR_threshold" = NA
                                      , "leading_indicator_type" = NA
                                      , "EWS_threshold" = NA
                                      , "earliest_tp_EWS_date" = NA
                                      , "lead_lag_days" = NA
                                      , "ews_false_positive_n" = NA
                                      )

#par(mfrow = c(24,6))
#par(mfrow = c(1,1))
#par(mar=c(5,5,5,5))
#    plot( 0,0 
#      , xlim = c(-25,75)
#      , ylim = c(0,20)
#      , xlab = "EWS lead(-ve)/lag(+ve) (days)"
#      , ylab = "Number of false positive EWS"
#    )
counter = 0

for ( wave_no in sort( unique( wave_results_df$wave_n ) ) ){
  for ( lead_ind in unique( wave_results_df$leading_indicator_type ) ){
    for ( min_age in unique( wave_results_df$TFPS_cluster_min_age ) ){
      for ( max_age in unique( wave_results_df$TFPS_cluster_max_age ) ){
        for ( min_desc in unique( wave_results_df$TFPS_cluster_min_descendants ) ){
          for ( p_val_th in p_val_thresholds ){
            for ( lgr_th in parent_sub_lgr_threshold ){
              #' Loop through EWS threshold levels, which represent the classifier boundaries
              for ( ews_th in unique( wave_results_df$EWS_threshold ) ){
                #message( wave_no," ",lead_ind," ", min_age," ",max_age," ",min_desc," ",p_val_th," ",ews_th," ",lgr_th )
                #for ( i  in 1 : length( ews_thresholds ) ){
                #' Filter by leading indicator, scan variables, cluster filter (LGR p-value) and EWS threshold
                wave_results_df_filtered = data.frame()
                wave_results_df_filtered = subset( wave_results_df , wave_results_df$wave_n == wave_no &
                                                                     wave_results_df$leading_indicator_type == lead_ind &
                                                                     wave_results_df$TFPS_cluster_min_age == min_age &
                                                                     wave_results_df$TFPS_cluster_max_age == max_age &
                                                                     wave_results_df$TFPS_cluster_min_descendants == min_desc &
                                                                     wave_results_df$LGR_p_value_threshold == p_val_th &
                                                                     wave_results_df$parent_sub_LGR_threshold == lgr_th &
                                                                     wave_results_df$EWS_threshold == ews_th
                                                  )
                
                #'  Compile new row for new dataframe
                if ( nrow( wave_results_df_filtered ) > 0 ){
                  new_row = wave_results_analysis_df[ 1 , ]
                  
                  new_row$wave_n = wave_no
                  new_row$wave_name = wave_results_df_filtered$wave_name[1]
                  new_row$pango = wave_results_df_filtered$pango[1]
                  new_row$wave_inflection_start_date = wave_results_df_filtered$wave_inflection_start_date[1]
                  new_row$wave_band_start = wave_results_df_filtered$wave_band_start[1]
                  new_row$wave_band_end = wave_results_df_filtered$wave_band_end[1]
                  new_row$leading_indicator_type = lead_ind
                  new_row$TFPS_cluster_min_age = min_age
                  new_row$TFPS_cluster_max_age = max_age
                  new_row$TFPS_cluster_min_descendants = min_desc
                  new_row$LGR_p_value_threshold = p_val_th
                  new_row$parent_sub_LGR_threshold = lgr_th
                  new_row$EWS_threshold = ews_th
                  #' Calculate number of false positive EWS for each wave
                  new_row$ews_false_positive_n = nrow( subset( wave_results_df_filtered , wave_results_df_filtered$positive_EWS == FALSE ) )
                  #' Determine earliest true positive for each wave
                  if( nrow( subset ( wave_results_df_filtered , wave_results_df_filtered$positive_EWS == TRUE ) ) > 0 ){
                    new_row$earliest_tp_EWS_date = min( subset ( wave_results_df_filtered , wave_results_df_filtered$positive_EWS == TRUE )$EWS_date )
                    new_row$lead_lag_days = as.Date( new_row$earliest_tp_EWS_date , origin = "1970-01-01" ) - new_row$wave_inflection_start_date
                  } else { 
                    new_row$earliest_tp_EWS_date = NA
                    new_row$lead_lag_days = NA
                  }
                  
                  wave_results_analysis_df = rbind( wave_results_analysis_df , new_row )
              }
            }
            #' Plot relationship between number of false positives and EWS lead/lag
            #df = subset( wave_results_analysis_df , wave_results_analysis_df$wave_n == wave_no &
            #               wave_results_analysis_df$leading_indicator_type == lead_ind &
            #               wave_results_analysis_df$TFPS_cluster_min_age == min_age &
            #               wave_results_analysis_df$TFPS_cluster_max_age == max_age &
            #               wave_results_analysis_df$TFPS_cluster_min_descendants == min_desc &
            #               wave_results_analysis_df$LGR_p_value_threshold == p_val_th &
            #               wave_results_analysis_df$parent_sub_lgr_threshold == lgr_th
            #)
            #if ( sum( !is.na( df$lead_lag_days ) ) > 0 ){
            #  df = data.table::setorder( df , "lead_lag_days" )
            #  counter=counter+1
            #  x = df$lead_lag_days
            #  y = df$ews_false_positive_n
            #  #lines( x , y , col = counter , pch = 16 )
            #  fit1 <- lm( y ~ x , data=df )
            #  x_axis = seq(-25,75,1)
            #  lines(x_axis, predict(fit1, data.frame(x=x_axis)), col=counter)
            }
          }
        }
      }
    }
  }
}

#View( wave_results_analysis_df )
#' save dataframe
#setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_11/analysis")
saveRDS( wave_results_analysis_df , file = paste0( "wave_results_analysis_df_lgr_th_",parent_sub_lgr_threshold,"_pval_th_",p_val_thresholds,"_li_",filenames_prefix[ as.numeric( i ) ],".rds" ) )
#**WILL BE TOO BIG TO SAVE A .CSV**
#write.csv( wave_results_analysis_df , file = "wave_results_analysis_df.csv")