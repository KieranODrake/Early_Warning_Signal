#' Code to process/analyse output from tfp scanner. The product will be a data 
#' frame containing a time series of the variance of the cluster logistic growth
#' rates as well as other information. A list of lists will also be produced 
#' giving the individual cluster growth rates for each date.
#' 
#' 1 - remove clusters without a recent sample/sequence
#' 2 - remove clusters with sub-clusters
#' 3 - remove clusters with overlapping tips
#' 4 - record logistic growth rates of remaining clusters 
#' 5 - calculate variance of remaining cluster growth rates
#' 6 - output summary time series in a data frame

#' List of directories to check for tfp scan files
#path_base <- ""
dir_list <- c(    "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_08_min_desc_20/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_08_min_desc_30/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_08_min_desc_40/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_08_min_desc_50/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_08_min_desc_60/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_08_min_desc_70/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_08_min_desc_80/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_08_min_desc_90/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_08_min_desc_100/HPC/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_12_min_desc_20/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_12_min_desc_30/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_12_min_desc_40/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_12_min_desc_50/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_12_min_desc_60/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_12_min_desc_70/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_12_min_desc_80/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_12_min_desc_90/scan_outputs"
                , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/max_age_12_min_desc_100/scan_outputs"
             )
dir_name <- c(     "ma08md20"
                 , "ma08md30"
                 , "ma08md40"
                 , "ma08md50"
                 , "ma08md60"
                 , "ma08md70"
                 , "ma08md80"
                 , "ma08md90"
                 , "ma08md100"
                 , "ma12md20"
                 , "ma12md30"
                 , "ma12md40"
                 , "ma12md50"
                 , "ma12md60"
                 , "ma12md70"
                 , "ma12md80"
                 , "ma12md90"
                 , "ma12md100"
            )

dir_df <- data.frame( "name" = dir_name , "path" = dir_list )
rm( dir_list , dir_name )
#' Cycle through list of directories that contain tfp scan output files for 
#' analysis
for (m in 1 : nrow( dir_df ) ){
  message( "Analysing ", dir_df$name[ m ] )
  #' Select directory containing tfp scanner output 
  #' (the .rds files but not the env files)
  #tfps_dir = choose.dir()
  setwd( dir_df$path[ m ] )
  tfps_output_list <- list.files()
  
  #' Initiate vectors which will be added to iteratively and then formed into a 
  #' data frame after the for loop below
  date_list <- c()
  grth_rate_var_list <- c()
  grth_rate_min_list <- c()
  grth_rate_max_list <- c()
  grth_rate_mean_list <- c()
  grth_rate_wtd_mean_list <- c()
  n_clusters_raw_list <- c()
  n_clusters_na_grth_rate_list <- c()
  n_clusters_too_old_list <- c()
  n_clusters_sub_list <- c()
  n_clusters_overlap_list <- c()
  n_clusters_remaining_list <- c()
  
  logistic_growth_rate_list <- list()
  logistic_growth_rate_dates <- list()
  
  #start_time_total = Sys.time()
  #' Loop through all tfp scans in folder to create time series
  for ( i in 1: length( tfps_output_list ) ) {
    #start_time_single = Sys.time()
    message(i)
    tfps_output_filename = tfps_output_list[ i ] #"scanner-2020-08-14.rds"
    tfps_output = readRDS( tfps_output_filename )
    n_cluster_raw <- nrow( tfps_output )
    tfps_date = as.Date( substr( tfps_output_filename , 9 , 20 ) )
    tfps_output <- subset( tfps_output , !is.na( tfps_output$logistic_growth_rate ))
    n_cluster_na_grth_rate <- n_cluster_raw - nrow( tfps_output )
    
    tfps_output_filtered <- tfps_output
    
    #' 1 - Remove clusters without a recent sample/sequence (14 days)
    #cut_off <-  max( tfps_output$most_recent_tip ) - 14
    
    #tfps_output_filtered <- subset( tfps_output_filtered ,
    #                                tfps_output_filtered$most_recent_tip >= cut_off
    #                              )
    #n_cluster_too_old <- n_cluster_raw - n_cluster_na_grth_rate - nrow( tfps_output_filtered )
    #'**OR (if want to switch off the filtering out clusters without recent samples)
    n_cluster_too_old <- 0
    
    #' 2 - Remove clusters containing sub-clusters (and the sub-clusters)
    sub_clust <- intersect( tfps_output_filtered$cluster_id 
                            , tfps_output_filtered$parent_number )
    '%ni%' <- Negate("%in%") #' Define 'not in'
    tfps_output_filtered <- subset( tfps_output_filtered 
                                    , ( tfps_output_filtered$cluster_id %ni% sub_clust ) &
                                      ( tfps_output_filtered$parent_number %ni% sub_clust ) )
    n_cluster_sub <- n_cluster_raw - n_cluster_na_grth_rate - n_cluster_too_old - 
                                                      nrow( tfps_output_filtered )
    
    #' 3 - Remove clusters with overlapping tips
    
    #' Check if there are any overlapping tips before trying to identify the 
    #' overlapping clusters (which is more time consuming)
    tips_all = data.frame( "tips" = c())
    
    for ( j in 1 : nrow( tfps_output_filtered ) ) { 
      tips_single = as.data.frame( strsplit( tfps_output_filtered$tips[ j ] , "|" , fixed = TRUE ) )
      colnames(tips_single) <- "tips"
      tips_all <- as.data.frame( rbind( tips_all , tips_single ))
      #message( j , " of ", nrow( tfps_output_filtered ), " clusters")
    }
    if ( nrow( unique( tips_all ) ) < nrow( tips_all ) ){
      #message( "There is cluster overlapping" )
      overlap_assessment <- TRUE
    } else if ( nrow( unique( tips_all ) ) == nrow( tips_all ) ) { 
      #message( "There is NO cluster overlapping" )
      overlap_assessment <- FALSE
    }
    
    #' Remove clusters that have overlapping tips (remove both clusters). 
    #' Such clusters need to be removed as we only want to include 
    #' non-overlapping clusters in our analysis.
    
    if ( overlap_assessment == TRUE ){
      overlap_rows <- c()
      tip_freq <- data.frame( table( tips_all ) )
      tips_overlap <- data.frame( subset( tip_freq , tip_freq$Freq > 1 )$tips )
      colnames( tips_overlap ) <- "tips"
      for ( k in 1 : nrow( tfps_output_filtered ) ) {
        tips_k <- as.data.frame( strsplit( tfps_output_filtered$tips[ k ] , "|" , fixed = TRUE ) )
        colnames( tips_k ) <- "tips"
        n_overlaps = nrow( dplyr::intersect( tips_overlap , tips_k ) )
        if ( n_overlaps > 0 ){
          overlap_rows <- c( overlap_rows , k )
        }
      }
      n_overlap_rows <- length( overlap_rows )
      tfps_output_filtered_overlap <- data.frame( tfps_output_filtered[ -overlap_rows , ] )
    } else if ( overlap_assessment == FALSE ) { 
      tfps_output_filtered_overlap = tfps_output_filtered
      n_overlap_rows = 0
      }
    
    #' 4 - Record set of logistical growth rates for clusters (after filtering) for each scan
    logistic_growth_rate_list[[ i ]] <- tfps_output_filtered_overlap$logistic_growth_rate
    
    #' 5 - Calculate variance of cluster logistical growth rates, weighted by cluster size
    number_clusters <- nrow( tfps_output_filtered_overlap )
    growth_rates <- tfps_output_filtered_overlap$logistic_growth_rate
    cluster_sizes <- tfps_output_filtered_overlap$cluster_size
    wtd_mean_growth <- sum( cluster_sizes * growth_rates ) / sum( cluster_sizes )
    cluster_growth_var <- sum( ( growth_rates - wtd_mean_growth ) ^ 2 ) / number_clusters
    
    #' Plots
    #hist(growth_rates,breaks = 30)
    #plot(density(growth_rates))
    #plot(growth_rates)
    
    #' Add scan analysis to variable lists
    date_list                    <- c( date_list , tfps_date )
    grth_rate_var_list           <- c( grth_rate_var_list , cluster_growth_var )
    grth_rate_min_list           <- c( grth_rate_min_list , min( growth_rates ) )
    grth_rate_max_list           <- c( grth_rate_max_list , max( growth_rates ) )
    grth_rate_mean_list          <- c( grth_rate_mean_list , mean( growth_rates) )
    grth_rate_wtd_mean_list      <- c( grth_rate_wtd_mean_list , wtd_mean_growth )
    n_clusters_raw_list          <- c( n_clusters_raw_list , n_cluster_raw )
    n_clusters_na_grth_rate_list <- c( n_clusters_na_grth_rate_list , n_cluster_na_grth_rate )
    n_clusters_too_old_list      <- c( n_clusters_too_old_list , n_cluster_too_old )
    n_clusters_sub_list          <- c( n_clusters_sub_list , n_cluster_sub )
    n_clusters_overlap_list      <- c( n_clusters_overlap_list , n_overlap_rows )
    n_clusters_remaining_list    <- c( n_clusters_remaining_list , number_clusters )
    
    #end_time_single = Sys.time()
    #print( i )
    #print( end_time_single - start_time_single )
  }
  #end_time_total = Sys.time()
  #print( i )
  #print( end_time_total - start_time_total )
  
  #' 6 - output summary data for time series of scans/trees into data frame
  #' **Possibly add quantiles to this data frame**
  tfps_growth_var <- data.frame(  "date"                    <- as.Date( date_list , origin = "1970-01-01")
                                , "growth_rate_var"         <- grth_rate_var_list
                                , "grth_rate_min"           <- grth_rate_min_list
                                , "grth_rate_max"           <- grth_rate_max_list
                                , "grth_rate_mean"          <- grth_rate_mean_list
                                , "grth_rate_wtd_mean"      <- grth_rate_wtd_mean_list
                                , "n_clusters_raw"          <- n_clusters_raw_list
                                , "n_clusters_na_grth_rate" <- n_clusters_na_grth_rate_list
                                , "n_clusters_too_old"      <- n_clusters_too_old_list
                                , "n_clusters_sub"          <- n_clusters_sub_list
                                , "n_clusters_overlap"      <- n_clusters_overlap_list
                                , "n_clusters_remaining"    <- n_clusters_remaining_list
                                )
  colnames( tfps_growth_var ) <- c("date"
                                   , "growth_rate_var"
                                   , "grth_rate_min"
                                   , "grth_rate_max"
                                   , "grth_rate_mean"
                                   , "grth_rate_wtd_mean"
                                   , "n_clusters_raw"
                                   , "n_clusters_na_grth_rate"
                                   , "n_clusters_too_old"
                                   , "n_clusters_sub"
                                   , "n_clusters_overlap"
                                   , "n_clusters_remaining"
                                  )
  
  names( logistic_growth_rate_list ) <- date_list # Note that these dates are in numeric format and will need to be converted using as.Date( date_list[i], origin = "1970-01-01")

  #' Rename tfps_growth_var data frame and logistic_growth_rate_list with relevant names
  new_df_name <- paste( "tfps_" , dir_df$name[ m ] , "_growth_var_df" , sep = "" )
  assign( new_df_name , tfps_growth_var )
  
  new_list_name <- paste( "tfps_" , dir_df$name[ m ] , "_growth_rate_list" , sep = "" )
  assign( new_list_name , logistic_growth_rate_list )
  #' Remove data frame and list with generic names so no confusion when create 
  #' for next data set
  rm( tfps_growth_var )
  rm( logistic_growth_rate_list )
  
}

#' Compile logistic growth rate variance for different scan variables into a 
#' single dataframe
tfps_lgr_df = data.frame(   "date"         = tfps_ma08md20_growth_var_df$date
                       , "ma_08_md_20"  = tfps_ma08md20_growth_var_df$growth_rate_var
                       , "ma_08_md_30"  = tfps_ma08md30_growth_var_df$growth_rate_var
                       , "ma_08_md_40"  = tfps_ma08md40_growth_var_df$growth_rate_var
                       , "ma_08_md_50"  = tfps_ma08md50_growth_var_df$growth_rate_var
                       , "ma_08_md_60"  = tfps_ma08md60_growth_var_df$growth_rate_var
                       , "ma_08_md_70"  = tfps_ma08md70_growth_var_df$growth_rate_var
                       #, "ma_08_md_80"  = tfps_ma08md80_growth_var_df$growth_rate_var
                       #, "ma_08_md_90"  = tfps_ma08md90_growth_var_df$growth_rate_var
                       , "ma_08_md_100" = tfps_ma08md100_growth_var_df$growth_rate_var
                       , "ma_12_md_20"  = tfps_ma12md20_growth_var_df$growth_rate_var
                       , "ma_12_md_30"  = tfps_ma12md30_growth_var_df$growth_rate_var
                       , "ma_12_md_40"  = tfps_ma12md40_growth_var_df$growth_rate_var
                       , "ma_12_md_50"  = tfps_ma12md50_growth_var_df$growth_rate_var
                       , "ma_12_md_60"  = tfps_ma12md60_growth_var_df$growth_rate_var
                       , "ma_12_md_70"  = tfps_ma12md70_growth_var_df$growth_rate_var
                       #, "ma_12_md_80"  = tfps_ma12md80_growth_var_df$growth_rate_var
                       #, "ma_12_md_90"  = tfps_ma12md90_growth_var_df$growth_rate_var
                       , "ma_12_md_100" = tfps_ma12md100_growth_var_df$growth_rate_var
                      )

folder = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/Analysis"
setwd( folder )
write.csv( tfps_lgr_df , file="tfps_lgr.csv")

tfps_growth_rate_lists = list(   ma08md20mca1 = tfps_ma08md20_growth_rate_list
                               , ma08md30mca1 = tfps_ma08md30_growth_rate_list
                               , ma08md40mca1 = tfps_ma08md40_growth_rate_list
                               , ma08md50mca1 = tfps_ma08md50_growth_rate_list
                               , ma08md60mca1 = tfps_ma08md60_growth_rate_list
                               , ma08md70mca1 = tfps_ma08md70_growth_rate_list
                               , ma08md80mca1 = tfps_ma08md80_growth_rate_list
                               , ma08md90mca1 = tfps_ma08md90_growth_rate_list
                               , ma08md100mca1 = tfps_ma08md100_growth_rate_list
                               , ma12md20mca1 = tfps_ma12md20_growth_rate_list
                               , ma12md30mca1 = tfps_ma12md30_growth_rate_list
                               , ma12md40mca1 = tfps_ma12md40_growth_rate_list
                               , ma12md50mca1 = tfps_ma12md50_growth_rate_list
                               , ma12md60mca1 = tfps_ma12md60_growth_rate_list
                               , ma12md70mca1 = tfps_ma12md70_growth_rate_list
                               , ma12md80mca1 = tfps_ma12md80_growth_rate_list
                               , ma12md90mca1 = tfps_ma12md90_growth_rate_list
                               , ma12md100mca1 = tfps_ma12md100_growth_rate_list
                              )  

folder = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_08/Analysis"
setwd( folder )
saveRDS( tfps_growth_rate_lists , "tfps_growth_lists.RData" )

###############################################################################

#' Charts
library("data.table")
# Load data, trim and organise
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/GitHub/Early_Warning_Signal'
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
dat_type <- "hospitalisations"
#' data_load function is in 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/GitHub/Early_Warning_Signal'
hosp_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)

#' Plot of variance of cluster growth rates
df_dates = data.frame(    tfps_ma08md20_growth_var_df$date
                        , tfps_ma08md30_growth_var_df$date
                        , tfps_ma08md40_growth_var_df$date
                        , tfps_ma08md50_growth_var_df$date
                        , tfps_ma08md60_growth_var_df$date
                        , tfps_ma08md70_growth_var_df$date
                        #, tfps_ma08md70_growth_var_df$date
                        #, tfps_ma08md90_growth_var_df$date
                        , tfps_ma08md100_growth_var_df$date
                        , tfps_ma12md20_growth_var_df$date
                        , tfps_ma12md30_growth_var_df$date
                        , tfps_ma12md40_growth_var_df$date
                        , tfps_ma12md50_growth_var_df$date
                        , tfps_ma12md60_growth_var_df$date
                        , tfps_ma12md70_growth_var_df$date
                        , tfps_ma12md70_growth_var_df$date
                        , tfps_ma12md90_growth_var_df$date
                        , tfps_ma12md100_growth_var_df$date
                        
                      )
df_lgr = data.frame(   tfps_ma08md20_growth_var_df$growth_rate_var
                       , tfps_ma08md30_growth_var_df$growth_rate_var
                       , tfps_ma08md40_growth_var_df$growth_rate_var
                       , tfps_ma08md50_growth_var_df$growth_rate_var
                       , tfps_ma08md60_growth_var_df$growth_rate_var
                       , tfps_ma08md70_growth_var_df$growth_rate_var
                       #, tfps_ma08md80_growth_var_df$growth_rate_var
                       #, tfps_ma08md90_growth_var_df$growth_rate_var
                       , tfps_ma08md100_growth_var_df$growth_rate_var
                       , tfps_ma12md20_growth_var_df$growth_rate_var
                       , tfps_ma12md30_growth_var_df$growth_rate_var
                       , tfps_ma12md40_growth_var_df$growth_rate_var
                       , tfps_ma12md50_growth_var_df$growth_rate_var
                       , tfps_ma12md60_growth_var_df$growth_rate_var
                       , tfps_ma12md70_growth_var_df$growth_rate_var
                       , tfps_ma12md80_growth_var_df$growth_rate_var
                       , tfps_ma12md90_growth_var_df$growth_rate_var
                       , tfps_ma12md100_growth_var_df$growth_rate_var
                     )

data_colour = c("green", "black", "red", "purple", "blue", "orange", "darkgreen"
                , "green", "black", "red", "purple", "blue", "orange", "darkgreen", "yellow")

plot( df_dates[ , 1 ]
      , df_lgr[ , 1 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      , xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 1.6 ) 
      , col = data_colour[ 1 ])

for ( i in 2 : length( df_dates ) ){
  lines( df_dates[ , i ]
         , df_lgr[ , i ]
         , typ="l"
         , col=data_colour[ i ])
}

points( hosp_df$date
        , hosp_df$cases/2000
        #, typ="l"
        , col="grey")

#legend("topright",legend=c("max age = 8wks, min descendants = 20"
#                           ,"max age = 8wks, min descendants = 100"
#                           ,"max age = 12wks, min descendants = 20"
#                           ,"max age = 12wks, min descendants = 100"
#                           ,"UK Covid-19 hospitalisations (1/5000)")
#       ,lty=c(1,1,1,1,1)
#       ,col=c("red","green","blue","orange","black"))

###############################
plot( lgr_ma_08_md_100$date
      , lgr_ma_08_md_100$growth_rate_var
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      , ylim = c( 0 , 1.6 ) 
      , col="green")

lines( lgr_ma_08_md_20$date
       , lgr_ma_08_md_20$growth_rate_var
       , typ="l"
       , col="red")
lines( lgr_ma_12_md_20$date
       , lgr_ma_12_md_20$growth_rate_var
       , typ="l"
       , col="blue")
lines( lgr_ma_12_md_100$date
       , lgr_ma_12_md_100$growth_rate_var
       , typ="l"
       , col="orange")
points( hosp_df$date
       , hosp_df$cases/5000
       #, typ="l"
       , col="black")

legend("topright",legend=c("max age = 8wks, min descendants = 20"
                           ,"max age = 8wks, min descendants = 100"
                           ,"max age = 12wks, min descendants = 20"
                           ,"max age = 12wks, min descendants = 100"
                           ,"UK Covid-19 hospitalisations (1/5000)")
       ,lty=c(1,1,1,1,1)
       ,col=c("red","green","blue","orange","black"))

#' Number of clusters after filtering
plot( lgr_ma_08_md_100$date
      , lgr_ma_08_md_100$n_clusters_remaining
      , typ = "l" 
      , xlab = "date"
      , ylab = "Number of clusters after filtering"
      , ylim = c( 0 , 400 ) 
      , col="green") 
lines( lgr_ma_08_md_20$date
       , lgr_ma_08_md_20$n_clusters_remaining
       , typ="l"
       , col="red")
lines( lgr_ma_12_md_20$date
       , lgr_ma_12_md_20$n_clusters_remaining
       , typ="l"
       , col="blue")
lines( lgr_ma_12_md_100$date
       , lgr_ma_12_md_100$n_clusters_remaining
       , typ="l"
       , col="orange")
points( hosp_df$date
        , hosp_df$cases/20
        #, typ="l"
        , col="black")

legend("topleft",legend=c("max age = 8wks, min descendants = 20"
                           ,"max age = 8wks, min descendants = 100"
                           ,"max age = 12wks, min descendants = 20"
                           ,"max age = 12wks, min descendants = 100"
                           ,"UK Covid-19 hospitalisations (1/20)")
       ,lty=c(1,1,1,1,1)
       ,col=c("red","green","blue","orange","black"))



##############################################################################
#' Comparison of number of clusters pre and post filtering
lgr = tfps_ma08md100_growth_var_df
par(mar = c(5, 4, 4, 4))
#' Plot hospitalisation data for context
plot(   x = hosp_df$date
        , y = hosp_df$cases
        , typ="l"
        , xlab = "date"
        , ylab = ""
        #, ylab = "number of clusters"
        , col="grey")
#' All clusters output from tfp scanner
lines(   x = lgr$date
      , y = lgr$n_clusters_remaining + lgr$n_clusters_overlap + lgr$n_clusters_sub + lgr$n_clusters_too_old + lgr$n_clusters_na_grth_rate
      #, y = lgr$n_clusters_raw
      , typ="l"
      , col = "black")
#' Filter out clusters with no calculated logistic growth rate
lines(  x = lgr$date
        , y = lgr$n_clusters_remaining + lgr$n_clusters_overlap + lgr$n_clusters_sub + lgr$n_clusters_too_old
        , typ ="l"
        , col = "orange" )
#' Filter out non-extant clusters i.e. those without a sequence within 2 weeks of the tree date
lines(  x = lgr$date
        , y = lgr$n_clusters_remaining + lgr$n_clusters_overlap + lgr$n_clusters_sub
        , typ ="l"
        , col = "dark green" )
#' Filter out sub-clusters (both parent and sub-cluster)
lines(  x = lgr$date
        , y = lgr$n_clusters_remaining + lgr$n_clusters_overlap
        , typ ="l"
        , col = "purple" )
#' Filter out clusters with overlapping (matching) tips (both clusters removed)
lines(  x = lgr$date
        , y = lgr$n_clusters_remaining
        , typ ="l"
        , col = "red" )
#' Add line at y = 0
lines(x = hosp_df$date
      , y = replicate( length( hosp_df$date ) , 0 )
      , typ ="l"
      , col = "black" )
#' Add legend
legend("topleft"
       , legend = c( "All clusters from tfp scan" 
                     , "Removed: no logistic growth rate"
                     , "Removed: non-extant (2 weeks)"
                     , "Removed: sub-clusters (and their parents)"
                     , "Removed: overlapping tips (both)"
                     , "Covid-19 hospitalisations"
                     )
       , lty = c(1,1,1,1,1,1)
       , col = c("black","orange","darkgreen","purple","red","grey")
       , cex = 0.75)

################
#' Similar plot but with percentages rather than absolute number of clusters

par(mar=c(5, 5, 3, 17), xpd=FALSE)
#' All clusters output from tfp scanner
plot(   x = lgr$date
         , y = 100 * ( ( lgr$n_clusters_remaining + lgr$n_clusters_overlap + lgr$n_clusters_sub + lgr$n_clusters_too_old + lgr$n_clusters_na_grth_rate ) / lgr$n_clusters_raw )
         #, y = lgr$n_clusters_raw
        , xlab = "date"
        , ylab = "% of total clusters"
        , ylim = c(0,100)
        , typ="l"
         , col = "black")
#' Filter out clusters with no calculated logistic growth rate
lines(  x = lgr$date
        , y = 100 * ( ( lgr$n_clusters_remaining + lgr$n_clusters_overlap + lgr$n_clusters_sub + lgr$n_clusters_too_old ) / lgr$n_clusters_raw )
        , typ ="l"
        , col = "orange" )
#' Filter out non-extant clusters i.e. those without a sequence within 2 weeks of the tree date
lines(  x = lgr$date
        , y = 100 * ( ( lgr$n_clusters_remaining + lgr$n_clusters_overlap + lgr$n_clusters_sub ) / lgr$n_clusters_raw )
        , typ ="l"
        , col = "dark green" )
#' Filter out sub-clusters (both parent and sub-cluster)
lines(  x = lgr$date
        , y = 100 * ( ( lgr$n_clusters_remaining + lgr$n_clusters_overlap) / lgr$n_clusters_raw )
        , typ ="l"
        , col = "purple" )
#' Filter out clusters with overlapping (matching) tips (both clusters removed)
lines(  x = lgr$date
        , y = 100 * (lgr$n_clusters_remaining / lgr$n_clusters_raw )
        , typ ="l"
        , col = "red" )
#' Add line at y = 0
lines(x = lgr$date
      , y = replicate( length( lgr$date ) , 0 )
      , typ ="l"
      , col = "black" )
lines(   x = hosp_df$date
        , y = hosp_df$cases / 50
        , typ="l"
        #, ylab = "number of clusters"
        , col="grey")
#' Add legend
par(mar=c(5, 5, 3, 17), xpd=TRUE)
legend("topright"
       , inset = c( -0.6, 0)
       , legend = c( "All clusters from tfp scan" 
                     , "Removed: no logistic growth rate"
                     , "Removed: non-extant (2 weeks)"
                     , "Removed: sub-clusters (and their parents)"
                     , "Removed: overlapping tips (both)"
                     , "Covid-19 hospitalisations (1/40)"
       )
       , lty = c(1,1,1,1,1,1)
       , col = c("black","orange","darkgreen","purple","red","grey")
       , cex = 0.75)

################################################################################
#'Comparing variance in growth rates for clusters including and excluding
#'those without a sample within 2 weeks of the tree date 

#' Charts
library("data.table")
# Load data, trim and organise
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/GitHub/Early_Warning_Signal'
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
dat_type <- "hospitalisations"
#' data_load function is in 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/GitHub/Early_Warning_Signal'
hosp_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)

#' Plot of variance of cluster growth rates
par(mar = c(5, 4, 4, 2) + 0.1 , xpd=FALSE )
par(mfrow = c(1,1))

plot( tfps_lgr_df[ , 1 ]
      , tfps_lgr_df[ , 2 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 1.6 ) 
      , col = "black")

lines( tfps_lgr_non_extant_df[ , 1 ]
       , tfps_lgr_non_extant_df[ , 2 ]
         , typ="l"
         , col="red")

lines( tfps_lgr_df[ , 1 ]
       , tfps_lgr_df[ , 9 ]
       , typ="l"
       , col="blue")

points( hosp_df$date
        , hosp_df$cases/3000
        #, typ="l"
        , col="grey")

legend("topright"
       #, inset = c( -0.6, 0)
       , legend = c( "All filters applied" 
                     , "No removal of non-extant clusters"
                     , "Covid-19 hospitalisations (1/3000)"
       )
       , lty = c(1,1,1)
       , col = c("black","red","grey")
       #, cex = 0.75
       )

#' plot difference
plot( tfps_lgr_df[ , 1 ]
      , tfps_lgr_non_extant_df[ , 2 ] - tfps_lgr_df[ , 2 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "Difference in variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      #, ylim = c( 0 , 1.6 ) 
      , col = "blue")
lines( tfps_lgr_df[ , 1 ] , replicate(nrow(tfps_lgr_df),0),col="black")
