#' Code to process/analyse output from tfp scanner. The product will be a data 
#' frame containing a time series of the variance of the cluster logistic growth
#' rates as well as other information. A list of lists will also be produced 
#' giving the individual cluster growth rates for each date.
#' 
#' 1 - remove clusters without a recent sample/sequence
#' 2 - remove clusters with sub-clusters
#' 3 - remove clusters with overlapping tips
#' 4 - remove clusters that are not external
#' 5.1 - record logistic growth rates of remaining clusters
#' 5.2 - record cluster sizes of remaining clusters 
#' 5.3 - record clock outlier statistic
#' 6 - calculate variance of remaining cluster growth rates
#' 7 - output summary time series in data frames

library( magrittr ) 
library( lubridate ) 
library( ggplot2 )
library( stringr )
library(zoo)

#' List of variables
min_age = c( 7 , 14 , 28 )
max_age = c( 56 , 84 )
min_desc = c( 20 , 50 , 100 )

#' Determine whether the minimum number of descendants (cluster size) is fixed 
#' or proportional to the number of samples in during the max-age period
#min_desc_type = "fixed"  
min_desc_type = "proportional"

#' Recreate file name suffix
counter = 0
fn_suffix = c()
for (i in 1 : length( min_age ) ){
  for (j in 1:length( max_age )){
    #' Next step depends on whether min descendants was fixed (and so consistently named in filenames)
    #' or is different depending on the date of the tree
    if ( min_desc_type == "fixed"){
      for (k in 1:length( min_desc )){
        counter = counter + 1
        fn_suffix[ counter ] = paste( "min_age_" , min_age[ i ] , "-max_age_" , max_age[ j ] , "-min_desc_" , min_desc[ k ], sep="" )
      }
    }
    if ( min_desc_type == "proportional" ){
      counter = counter + 1
      fn_suffix[ counter ] = paste( "min_age_" , min_age[ i ] , "-max_age_" , max_age[ j ] , sep="" )
    }
  }
}


rm( counter , i , j , k )

#' TFPScan file output of format 'scanner-2021-12-22-min_age_7-max_age_84-min_desc_20.rds'
#' 
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/outputs_perc/")
#setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/outputs/")
#fns = list.files( 'outputs', full.name=TRUE )
fns = list.files()

#' Select filters to apply to clusters (TRUE in each case means that filter will be applied)
#'**STILL NEED TO BE INCORPORATED INTO MAIN CODE**
#'#' 1 = remove clusters without a recent sample/sequence
extant = TRUE
#' 2 = remove clusters with sub-clusters
non_sub_clusters = TRUE
#' 3 = remove clusters with overlapping tips
non_overlap = TRUE
#' 4 = remove clusters that are not external
external = TRUE

#' Cycle through sets of variables for analysis
for ( n in 1 : length( fn_suffix ) ) {
  if ( min_desc_type == "fixed" ){ file_str_end = 51 }
  if ( min_desc_type == "proportional" ){ file_str_end = 39 }
  #' This matches the middle part of the TFPScan filename 'scanner-2021-12-22-min_age_7-max_age_84-min_desc_20.rds'  
  file_list = subset( fns , substr( fns , 20 , file_str_end ) == fn_suffix[ n ] )
  #' if file_list returns nothing then change the length of the character string to match
  if ( length( file_list ) == 0 ){
    file_list = subset( fns , substr( fns , 20 , file_str_end + 1 ) == fn_suffix[ n ] )
  }
  if ( length( file_list ) == 0 ){
    file_list = subset( fns , substr( fns , 20 , file_str_end + 2 ) == fn_suffix[ n ] )
  }


  #' Create list of scan dates
  scandates <- as.Date( regmatches( file_list , regexpr('(\\d{4}-\\d{2}-\\d{2})', file_list) ) )

  #' (Re)initiate vectors which will be added to iteratively and then formed into a 
  #' data frame after the for loop below
  date_list <- c()
  #'
  grth_rate_var_wtd_list <- c()
  #' variance: population and sample for logistic growth rates, simple logistic growth rates and GAM logistic growth rates
  grth_rate_var_pop_list <- c()
  grth_rate_var_samp_list <- c()
  grth_rate_simple_var_pop_list <- c()
  grth_rate_simple_var_samp_list <- c()
  grth_rate_gam_var_pop_list <- c()
  grth_rate_gam_var_samp_list <- c()
  #' Minimum for logistic growth rates, simple logistic growth rates and GAM logistic growth rates
  grth_rate_min_list <- c()
  grth_rate_simple_min_list <- c()
  grth_rate_gam_min_list <- c()
  #' Maximum for logistic growth rates, simple logistic growth rates and GAM logistic growth rates
  grth_rate_max_list <- c()
  grth_rate_simple_max_list <- c()
  grth_rate_gam_max_list <- c()
  #' Mean for logistic growth rates, simple logistic growth rates and GAM logistic growth rates
  grth_rate_mean_list <- c()
  grth_rate_simple_mean_list <- c()
  grth_rate_gam_mean_list <- c()
  #' Mean Weighted by cluster size for logistic growth rates, simple logistic growth rates and GAM logistic growth rates
  grth_rate_wtd_mean_list <- c()
  grth_rate_simple_wtd_mean_list <- c()
  grth_rate_gam_wtd_mean_list <- c()
  #' clock outlier stats - calculated on absolute values as sign is not important
  clock_outlier_max_list <- c( )
  clock_outlier_mean_list <- c()
  #' Number of clusters remaining after various filters
  n_clusters_raw_list <- c()
  n_clusters_na_grth_rate_list <- c()
  n_clusters_too_old_list <- c()
  n_clusters_sub_list <- c()
  n_clusters_overlap_list <- c()
  n_clusters_remaining_list <- c()
  mean_cluster_size_list <- c()
  
  
  #n_clusters_above_min_desc_prop_56 <- c()
  #n_clusters_above_min_desc_prop_84 <- c()
  
  logistic_growth_rate_list        <- list()
  logistic_growth_rate_simple_list <- list()
  logistic_growth_rate_gam_list    <- list()
  logistic_growth_rate_max_list    <- list()
  cluster_size_list                <- list()
  clock_outlier_list               <- list()
  
    #start_time_total = Sys.time()
    #' Loop through all tfp scans in list to create time series
    for ( i in 1: length( file_list ) ) {
      #start_time_single = Sys.time()
      message("Analysing dataset " , n , ", file " , i , ": ", file_list[ i ] )
      tfps_output_filename = file_list[ i ] #"scanner-2020-08-14-min_age_7-max_age_56-min_desc_100.rds"
      tfps_output = readRDS( tfps_output_filename )
      n_cluster_raw <- nrow( tfps_output )
      tfps_date = as.Date( substr( tfps_output_filename , 9 , 20 ) )
      tfps_output <- subset( tfps_output , !is.na( tfps_output$logistic_growth_rate )) #$logistic_growth_rate or gam_logistic_growth_rate
      n_cluster_na_grth_rate <- n_cluster_raw - nrow( tfps_output )
    
      tfps_output_filtered <- tfps_output
    
      #' 1 - Remove clusters without a recent sample/sequence (14 days)
      cut_off <-  max( tfps_output$most_recent_tip ) - 14
      
      tfps_output_filtered <- subset( tfps_output_filtered ,
                                      tfps_output_filtered$most_recent_tip >= cut_off
                                    )
      n_cluster_too_old <- n_cluster_raw - n_cluster_na_grth_rate - nrow( tfps_output_filtered )
      #'**OR (if want to switch off the filtering out clusters without recent samples)
      #n_cluster_too_old <- 0
      
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
      if ( nrow( tfps_output_filtered ) > 0 ){ #' sometimes all clusters are removed due to sub-clusters - so need to check
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
      } else if ( nrow( tfps_output_filtered ) == 0 ){
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
      
      #' 4 - Remove non-external clusters from dataset
      tfps_output_filtered_overlap = subset( tfps_output_filtered_overlap , external_cluster == TRUE )
      
      #' 5.1 - Record set of logistic growth rates for clusters (after filtering) for each scan
      logistic_growth_rate_simple_list[[ i ]] <- tfps_output_filtered_overlap$simple_logistic_growth_rate
      logistic_growth_rate_gam_list[[ i ]] <- tfps_output_filtered_overlap$gam_logistic_growth_rate
      logistic_growth_rate_list[[ i ]] <- tfps_output_filtered_overlap$logistic_growth_rate
      
      #' 5.2 - Record set of cluster sizes for clusters (after filtering) for each scan
      cluster_size_list[[ i ]] <- tfps_output_filtered_overlap$cluster_size
      
      #' 5.3 - Record set of clock outlier statistics (after filtering) for each scan
      clock_outlier_list[[ i ]] <- tfps_output_filtered_overlap$clock_outlier
      
      #' 6 - Calculate variance of cluster logistic growth rates, weighted by cluster size
      #number_clusters <- nrow( tfps_output_filtered_overlap )
      cluster_sizes <- tfps_output_filtered_overlap$cluster_size
      #' Calculate for different logistic growth rates
      #growth_rates <- tfps_output_filtered_overlap$logistic_growth_rate
      growth_cluster_size_df = data.frame( "growth_rates" = tfps_output_filtered_overlap$logistic_growth_rate
                                           , "cluster_sizes" = cluster_sizes )
      growth_cluster_size_df = subset( growth_cluster_size_df , !is.na( growth_cluster_size_df$growth_rates ) )
      number_clusters = nrow( growth_cluster_size_df )
      growth_rates = growth_cluster_size_df$growth_rates
      cluster_sizes = growth_cluster_size_df$cluster_sizes
      wtd_mean_growth <- sum( cluster_sizes * growth_rates ) / sum( cluster_sizes )
      cluster_growth_var_pop <- sum( ( growth_rates - wtd_mean_growth ) ^ 2 ) / number_clusters
      cluster_growth_var_samp <- sum( ( growth_rates - wtd_mean_growth ) ^ 2 ) / ( number_clusters - 1 )
      #' simple logistic growth rates
      growth_cluster_size_simple_df = data.frame( "growth_rates" = tfps_output_filtered_overlap$simple_logistic_growth_rate
                                           , "cluster_sizes" = cluster_sizes )
      growth_cluster_size_simple_df = subset( growth_cluster_size_simple_df , !is.na( growth_cluster_size_simple_df$growth_rates ) )
      number_clusters_simple = nrow( growth_cluster_size_simple_df )
      growth_rates_simple = growth_cluster_size_simple_df$growth_rates
      cluster_sizes_simple = growth_cluster_size_df$cluster_sizes
      #growth_rates_simple <- tfps_output_filtered_overlap$simple_logistic_growth_rate #logistic_growth_rate or gam_logistic_growth_rate
      wtd_mean_growth_simple <- sum( cluster_sizes_simple * growth_rates_simple ) / sum( cluster_sizes_simple )
      cluster_growth_simple_var_pop  <- sum( ( growth_rates_simple - wtd_mean_growth_simple ) ^ 2 ) / number_clusters_simple
      cluster_growth_simple_var_samp <- sum( ( growth_rates_simple - wtd_mean_growth_simple ) ^ 2 ) / ( number_clusters_simple - 1 )
      #' GAM logistic growth rates
      #growth_rates_gam <- tfps_output_filtered_overlap$gam_logistic_growth_rate #logistic_growth_rate or gam_logistic_growth_rate
      growth_cluster_size_gam_df = data.frame( "growth_rates" = tfps_output_filtered_overlap$gam_logistic_growth_rate
                                                  , "cluster_sizes" = cluster_sizes )
      growth_cluster_size_gam_df = subset( growth_cluster_size_gam_df , !is.na( growth_cluster_size_gam_df$growth_rates ) )
      number_clusters_gam = nrow( growth_cluster_size_gam_df )
      growth_rates_gam = growth_cluster_size_gam_df$growth_rates
      cluster_sizes_gam = growth_cluster_size_gam_df$cluster_sizes
      
      wtd_mean_growth_gam <- sum( na.omit( cluster_sizes_gam * growth_rates_gam) ) / sum( cluster_sizes_gam )
      cluster_growth_gam_var_pop  <- sum( ( growth_rates_gam - wtd_mean_growth_gam ) ^ 2 ) / number_clusters_gam
      cluster_growth_gam_var_samp <- sum( ( growth_rates_gam - wtd_mean_growth_gam ) ^ 2 ) / ( number_clusters_gam - 1 )
      
      #' Weight by inverse of mean number of clusters to reduce effect of when 
      #' mean cluster size is smaller the variance of growth rates will be noisier
      cluster_growth_var_wtd <- cluster_growth_var_samp / mean( tfps_output_filtered_overlap$cluster_size , na.rm = TRUE)
      
      #' Record clock outlier stats - calculated on absolute values as the sign is not important
      clock_outlier_max   <- max(  abs( tfps_output_filtered_overlap$clock_outlier ) , na.rm = TRUE)
      clock_outlier_mean  <- mean( abs( tfps_output_filtered_overlap$clock_outlier ) , na.rm = TRUE)
      
      #' Plots
      #hist(growth_rates,breaks = 30)
      #plot(density(growth_rates))
      #plot(growth_rates)
      
      #' Add scan analysis to variable lists
      date_list                      <- c( date_list , tfps_date )
      grth_rate_var_wtd_list         <- c( grth_rate_var_wtd_list , cluster_growth_var_wtd )
      grth_rate_var_pop_list         <- c( grth_rate_var_pop_list , cluster_growth_var_pop )
      grth_rate_var_samp_list        <- c( grth_rate_var_samp_list , cluster_growth_var_samp )
      grth_rate_simple_var_pop_list  <- c( grth_rate_simple_var_pop_list , cluster_growth_simple_var_pop )
      grth_rate_simple_var_samp_list <- c( grth_rate_simple_var_samp_list , cluster_growth_simple_var_samp )
      
      grth_rate_gam_var_pop_list     <- c( grth_rate_gam_var_pop_list , cluster_growth_gam_var_pop )
      
      grth_rate_gam_var_samp_list    <- c( grth_rate_gam_var_samp_list , cluster_growth_gam_var_samp )
      
      grth_rate_min_list             <- c( grth_rate_min_list , min( growth_rates ) )
      grth_rate_simple_min_list      <- c( grth_rate_simple_min_list , min( growth_rates_simple ) )
      
      grth_rate_gam_min_list         <- c( grth_rate_gam_min_list , min( growth_rates_gam ) )
      
      grth_rate_max_list             <- c( grth_rate_max_list , max( growth_rates ) )
      grth_rate_simple_max_list      <- c( grth_rate_simple_max_list , max( growth_rates_simple ) )
      
      grth_rate_gam_max_list         <- c( grth_rate_gam_max_list , max( growth_rates_gam ) )
      
      grth_rate_mean_list            <- c( grth_rate_mean_list , mean( growth_rates ) )
      grth_rate_simple_mean_list     <- c( grth_rate_simple_mean_list , mean( growth_rates_simple ) )
      
      grth_rate_gam_mean_list        <- c( grth_rate_gam_mean_list , mean( growth_rates_gam ) )
      
      grth_rate_wtd_mean_list        <- c( grth_rate_wtd_mean_list , wtd_mean_growth )
      grth_rate_simple_wtd_mean_list <- c( grth_rate_simple_wtd_mean_list , wtd_mean_growth_simple )
      
      grth_rate_gam_wtd_mean_list    <- c( grth_rate_gam_wtd_mean_list , wtd_mean_growth_gam )
      
      clock_outlier_max_list         <- c( clock_outlier_max_list , clock_outlier_max )
      clock_outlier_mean_list        <- c( clock_outlier_mean_list , clock_outlier_mean )
      n_clusters_raw_list            <- c( n_clusters_raw_list , n_cluster_raw )
      n_clusters_na_grth_rate_list   <- c( n_clusters_na_grth_rate_list , n_cluster_na_grth_rate )
      n_clusters_too_old_list        <- c( n_clusters_too_old_list , n_cluster_too_old )
      n_clusters_sub_list            <- c( n_clusters_sub_list , n_cluster_sub )
      n_clusters_overlap_list        <- c( n_clusters_overlap_list , n_overlap_rows )
      n_clusters_remaining_list      <- c( n_clusters_remaining_list , number_clusters )
      mean_cluster_size_list         <- c( mean_cluster_size_list , mean( tfps_output_filtered_overlap$cluster_size ) )
      
      
      #end_time_single = Sys.time()
      #print( i )
      #print( end_time_single - start_time_single )
    }
    #end_time_total = Sys.time()
    #print( i )
    #print( end_time_total - start_time_total )
  
    #' 7 - output summary data for time series of scans/trees into data frame
    #' **Possibly add quantiles to this data frame**
    tfps_growth_var <- data.frame(  "date"                        <- as.Date( date_list , origin = "1970-01-01")
                                  , "growth_rate_var_wtd"         <- grth_rate_var_wtd_list
                                  , "growth_rate_var_pop"         <- grth_rate_var_pop_list
                                  , "growth_rate_var_samp"        <- grth_rate_var_samp_list
                                  , "growth_rate_simple_var_pop"  <- grth_rate_simple_var_pop_list
                                  , "growth_rate_simple_var_samp" <- grth_rate_simple_var_samp_list
                                  , "growth_rate_gam_var_pop"     <- grth_rate_gam_var_pop_list
                                  , "growth_rate_gam_var_samp"    <- grth_rate_gam_var_samp_list
                                  , "grth_rate_min"               <- grth_rate_min_list
                                  , "grth_rate_simple_min"        <- grth_rate_simple_min_list
                                  , "grth_rate_gam_min"           <- grth_rate_gam_min_list
                                  , "grth_rate_max"               <- grth_rate_max_list
                                  , "grth_rate_simple_max"        <- grth_rate_simple_max_list
                                  , "grth_rate_gam_max"           <- grth_rate_gam_max_list
                                  , "grth_rate_mean"              <- grth_rate_mean_list
                                  , "grth_rate_simple_mean"       <- grth_rate_simple_mean_list
                                  , "grth_rate_gam_mean"          <- grth_rate_gam_mean_list
                                  , "grth_rate_wtd_mean"          <- grth_rate_wtd_mean_list
                                  , "grth_rate_simple_wtd_mean"   <- grth_rate_simple_wtd_mean_list
                                  , "grth_rate_gam_wtd_mean"      <- grth_rate_gam_wtd_mean_list
                                  , "clock_outlier_max"           <- clock_outlier_max_list
                                  , "clock_outlier_mean"          <- clock_outlier_mean_list
                                  , "n_clusters_raw"              <- n_clusters_raw_list
                                  , "n_clusters_na_grth_rate"     <- n_clusters_na_grth_rate_list
                                  , "n_clusters_too_old"          <- n_clusters_too_old_list
                                  , "n_clusters_sub"              <- n_clusters_sub_list
                                  , "n_clusters_overlap"          <- n_clusters_overlap_list
                                  , "n_clusters_remaining"        <- n_clusters_remaining_list
                                  , "mean_cluster_size"           <- mean_cluster_size_list
                                  )
    colnames( tfps_growth_var ) <- c("date"
                                     , "growth_rate_var_wtd"         
                                     , "growth_rate_var_pop"        
                                     , "growth_rate_var_samp"        
                                     , "growth_rate_simple_var_pop"  
                                     , "growth_rate_simple_var_samp"  
                                     , "growth_rate_gam_var_pop"
                                     , "growth_rate_gam_var_samp"    
                                     , "grth_rate_min"               
                                     , "grth_rate_simple_min"        
                                     , "grth_rate_gam_min"           
                                     , "grth_rate_max"               
                                     , "grth_rate_simple_max"       
                                     , "grth_rate_gam_max"           
                                     , "grth_rate_mean"             
                                     , "grth_rate_simple_mean"       
                                     , "grth_rate_gam_mean"         
                                     , "grth_rate_wtd_mean"          
                                     , "grth_rate_simple_wtd_mean"   
                                     , "grth_rate_gam_wtd_mean"
                                     , "clock_outlier_max"
                                     , "clock_outlier_mean"
                                     , "n_clusters_raw"              
                                     , "n_clusters_na_grth_rate"     
                                     , "n_clusters_too_old"          
                                     , "n_clusters_sub"
                                     , "n_clusters_overlap"
                                     , "n_clusters_remaining"
                                     , "mean_cluster_size"
                                    )
    
    names( logistic_growth_rate_list ) <- date_list # Note that these dates are in numeric format and will need to be converted using as.Date( date_list[i], origin = "1970-01-01")
    names( logistic_growth_rate_simple_list ) <- date_list # Note that these dates are in numeric format and will need to be converted using as.Date( date_list[i], origin = "1970-01-01")
    names( logistic_growth_rate_gam_list ) <- date_list # Note that these dates are in numeric format and will need to be converted using as.Date( date_list[i], origin = "1970-01-01")
    names( cluster_size_list ) <- date_list # Note that these dates are in numeric format and will need to be converted using as.Date( date_list[i], origin = "1970-01-01")
    names( clock_outlier_list ) <- date_list # Note that these dates are in numeric format and will need to be converted using as.Date( date_list[i], origin = "1970-01-01")
    
    #' Rename tfps_growth_var data frame and logistic_growth_rate_list with relevant names
    var_values = gsub( "-" , "_" , fn_suffix[ n ] )
    new_df_name <- paste( "tfps_" , var_values , "_growth_var_df" , sep = "" )
    assign( new_df_name , tfps_growth_var )
    
    new_list_name <- paste( "tfps_" , var_values , "_growth_rate_list" , sep = "" )
    assign( new_list_name , logistic_growth_rate_list )

    new_list_name_simple <- paste( "tfps_" , var_values , "_growth_rate_simple_list" , sep = "" )
    assign( new_list_name_simple , logistic_growth_rate_simple_list )
    
    new_list_name_gam <- paste( "tfps_" , var_values , "_growth_rate_gam_list" , sep = "" )
    assign( new_list_name_gam , logistic_growth_rate_gam_list )

    new_list_name_co <- paste( "tfps_" , var_values , "_clock_outlier_list" , sep = "" )
    assign( new_list_name_co , clock_outlier_list )
    
    new_list_name_cs <- paste( "tfps_" , var_values , "_cluster_size_list" , sep = "" )
    assign( new_list_name_cs , cluster_size_list )
    
    #' Remove data frame and list with generic names so no confusion when create 
    #' for next data set
    rm( tfps_growth_var )
    rm( logistic_growth_rate_list ) ; rm( new_list_name )
    rm( logistic_growth_rate_simple_list ) ; rm( new_list_name_simple )
    rm( logistic_growth_rate_gam_list ) ; rm( new_list_name_gam )
    rm( clock_outlier_list ) ; rm( new_list_name_oc )
    rm( cluster_size_list ) ; rm( new_list_name_cs )
    
  }


#' Compile logistic growth rate variance for different scan variables into a 
#' single dataframe
#' tfps_lgr_mean_df or tfps_gam_lgr_df or tfps_lgr_max_df or tfps_lgr_wtd_mean_df 
#' or tfps_vlgr_s_df (sample var instead of population var), tfps_vlgr_p_df (population instead of sample variance)
tfps_vlgr_samp_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                         , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$growth_rate_var_samp
                         , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$growth_rate_var_samp
                         , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$growth_rate_var_samp
                         , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$growth_rate_var_samp
                         , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$growth_rate_var_samp
                         , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$growth_rate_var_samp
                         , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$growth_rate_var_samp
                         , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$growth_rate_var_samp
                         , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$growth_rate_var_samp
                         , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$growth_rate_var_samp
                         , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$growth_rate_var_samp
                         , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$growth_rate_var_samp
                         , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$growth_rate_var_samp
                         , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$growth_rate_var_samp
                         , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$growth_rate_var_samp
                         , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$growth_rate_var_samp
                         , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$growth_rate_var_samp
                         , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$growth_rate_var_samp
                         , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$growth_rate_var_samp
                         , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$growth_rate_var_samp
                         , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$growth_rate_var_samp
                         , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$growth_rate_var_samp
                         , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$growth_rate_var_samp
                         , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$growth_rate_var_samp
                         )
tfps_vlgr_simple_samp_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                 , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$growth_rate_simple_var_samp
                                 , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$growth_rate_simple_var_samp
                                 , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$growth_rate_simple_var_samp
                                 , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$growth_rate_simple_var_samp
                                 , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$growth_rate_simple_var_samp
                                 , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$growth_rate_simple_var_samp
                                 , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$growth_rate_simple_var_samp
                                 , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$growth_rate_simple_var_samp
                                 , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$growth_rate_simple_var_samp
                                 , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$growth_rate_simple_var_samp
                                 , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$growth_rate_simple_var_samp
                                 , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$growth_rate_simple_var_samp
                                 , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$growth_rate_simple_var_samp
                                 , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$growth_rate_simple_var_samp
                                 , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$growth_rate_simple_var_samp
                                 , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$growth_rate_simple_var_samp
                                 , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$growth_rate_simple_var_samp
                                 , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$growth_rate_simple_var_samp
                                 , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$growth_rate_simple_var_samp
                                 , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$growth_rate_simple_var_samp
                                 , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$growth_rate_simple_var_samp
                                 , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$growth_rate_simple_var_samp
                                 , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$growth_rate_simple_var_samp
                                 , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$growth_rate_simple_var_samp
)
tfps_vlgr_gam_samp_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                 , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$growth_rate_gam_var_samp
                                 , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$growth_rate_gam_var_samp
                                 , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$growth_rate_gam_var_samp
                                 , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$growth_rate_gam_var_samp
                                 , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$growth_rate_gam_var_samp
                                 , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$growth_rate_gam_var_samp
                                 , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$growth_rate_gam_var_samp
                                 , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$growth_rate_gam_var_samp
                                 , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$growth_rate_gam_var_samp
                                 , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$growth_rate_gam_var_samp
                                 , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$growth_rate_gam_var_samp
                                 , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$growth_rate_gam_var_samp
                                 , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$growth_rate_gam_var_samp
                                 , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$growth_rate_gam_var_samp
                                 , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$growth_rate_gam_var_samp
                                 , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$growth_rate_gam_var_samp
                                 , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$growth_rate_gam_var_samp
                                 , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$growth_rate_gam_var_samp
                                 , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$growth_rate_gam_var_samp
                                 , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$growth_rate_gam_var_samp
                                 , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$growth_rate_gam_var_samp
                                 , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$growth_rate_gam_var_samp
                                 , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$growth_rate_gam_var_samp
                                 , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$growth_rate_gam_var_samp
                                )
tfps_vlgr_pop_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                              , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$growth_rate_var_pop
                              , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$growth_rate_var_pop
                              , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$growth_rate_var_pop
                              , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$growth_rate_var_pop
                              , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$growth_rate_var_pop
                              , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$growth_rate_var_pop
                              , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$growth_rate_var_pop
                              , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$growth_rate_var_pop
                              , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$growth_rate_var_pop
                              , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$growth_rate_var_pop
                              , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$growth_rate_var_pop
                              , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$growth_rate_var_pop
                              , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$growth_rate_var_pop
                              , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$growth_rate_var_pop
                              , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$growth_rate_var_pop
                              , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$growth_rate_var_pop
                              , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$growth_rate_var_pop
                              , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$growth_rate_var_pop
                              , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$growth_rate_var_pop
                              , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$growth_rate_var_pop
                              , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$growth_rate_var_pop
                              , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$growth_rate_var_pop
                              , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$growth_rate_var_pop
                              , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$growth_rate_var_pop
)
tfps_vlgr_simple_pop_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$growth_rate_simple_var_pop
                                , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$growth_rate_simple_var_pop
                                , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$growth_rate_simple_var_pop
                                , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$growth_rate_simple_var_pop
                                , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$growth_rate_simple_var_pop
                                , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$growth_rate_simple_var_pop
                                , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$growth_rate_simple_var_pop
                                , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$growth_rate_simple_var_pop
                                , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$growth_rate_simple_var_pop
                                , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$growth_rate_simple_var_pop
                                , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$growth_rate_simple_var_pop
                                , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$growth_rate_simple_var_pop
                                , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$growth_rate_simple_var_pop
                                , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$growth_rate_simple_var_pop
                                , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$growth_rate_simple_var_pop
                                , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$growth_rate_simple_var_pop
                                , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$growth_rate_simple_var_pop
                                , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$growth_rate_simple_var_pop
                                , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$growth_rate_simple_var_pop
                                , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$growth_rate_simple_var_pop
                                , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$growth_rate_simple_var_pop
                                , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$growth_rate_simple_var_pop
                                , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$growth_rate_simple_var_pop
                                , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$growth_rate_simple_var_pop
)
tfps_vlgr_gam_pop_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$growth_rate_gam_var_pop
                                , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$growth_rate_gam_var_pop
                                , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$growth_rate_gam_var_pop
                                , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$growth_rate_gam_var_pop
                                , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$growth_rate_gam_var_pop
                                , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$growth_rate_gam_var_pop
                                , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$growth_rate_gam_var_pop
                                , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$growth_rate_gam_var_pop
                                , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$growth_rate_gam_var_pop
                                , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$growth_rate_gam_var_pop
                                , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$growth_rate_gam_var_pop
                                , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$growth_rate_gam_var_pop
                                , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$growth_rate_gam_var_pop
                                , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$growth_rate_gam_var_pop
                                , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$growth_rate_gam_var_pop
                                , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$growth_rate_gam_var_pop
                                , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$growth_rate_gam_var_pop
                                , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$growth_rate_gam_var_pop
                                , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$growth_rate_gam_var_pop
                                , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$growth_rate_gam_var_pop
                                , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$growth_rate_gam_var_pop
                                , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$growth_rate_gam_var_pop
                                , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$growth_rate_gam_var_pop
                                , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$growth_rate_gam_var_pop
)
tfps_vlgr_wtd_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                              , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$growth_rate_var_wtd
                              , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$growth_rate_var_wtd
                              , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$growth_rate_var_wtd
                              , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$growth_rate_var_wtd
                              , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$growth_rate_var_wtd
                              , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$growth_rate_var_wtd
                              , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$growth_rate_var_wtd
                              , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$growth_rate_var_wtd
                              , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$growth_rate_var_wtd
                              , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$growth_rate_var_wtd
                              , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$growth_rate_var_wtd
                              , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$growth_rate_var_wtd
                              , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$growth_rate_var_wtd
                              , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$growth_rate_var_wtd
                              , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$growth_rate_var_wtd
                              , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$growth_rate_var_wtd
                              , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$growth_rate_var_wtd
                              , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$growth_rate_var_wtd
                              , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$growth_rate_var_wtd
                              , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$growth_rate_var_wtd
                              , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$growth_rate_var_wtd
                              , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$growth_rate_var_wtd
                              , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$growth_rate_var_wtd
                              , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$growth_rate_var_wtd
)
tfps_lgr_max_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                               , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$grth_rate_max
                               , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$grth_rate_max
                               , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$grth_rate_max
                               , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$grth_rate_max
                               , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$grth_rate_max
                               , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$grth_rate_max
                               , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$grth_rate_max
                               , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$grth_rate_max
                               , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$grth_rate_max
                               , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$grth_rate_max
                               , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$grth_rate_max
                               , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$grth_rate_max
                               , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$grth_rate_max
                               , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$grth_rate_max
                               , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$grth_rate_max
                               , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$grth_rate_max
                               , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$grth_rate_max
                               , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$grth_rate_max
                               , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$grth_rate_max
                               , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$grth_rate_max
                               , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$grth_rate_max
                               , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$grth_rate_max
                               , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$grth_rate_max
                               , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$grth_rate_max
)
tfps_lgr_simple_max_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                               , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$grth_rate_simple_max
                               , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$grth_rate_simple_max
                               , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$grth_rate_simple_max
                               , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$grth_rate_simple_max
                               , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$grth_rate_simple_max
                               , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$grth_rate_simple_max
                               , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$grth_rate_simple_max
                               , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$grth_rate_simple_max
                               , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$grth_rate_simple_max
                               , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$grth_rate_simple_max
                               , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$grth_rate_simple_max
                               , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$grth_rate_simple_max
                               , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$grth_rate_simple_max
                               , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$grth_rate_simple_max
                               , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$grth_rate_simple_max
                               , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$grth_rate_simple_max
                               , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$grth_rate_simple_max
                               , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$grth_rate_simple_max
                               , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$grth_rate_simple_max
                               , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$grth_rate_simple_max
                               , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$grth_rate_simple_max
                               , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$grth_rate_simple_max
                               , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$grth_rate_simple_max
                               , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$grth_rate_simple_max
)
tfps_lgr_gam_max_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                      , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$grth_rate_gam_max
                                      , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$grth_rate_gam_max
                                      , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$grth_rate_gam_max
                                      , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$grth_rate_gam_max
                                      , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$grth_rate_gam_max
                                      , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$grth_rate_gam_max
                                      , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$grth_rate_gam_max
                                      , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$grth_rate_gam_max
                                      , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$grth_rate_gam_max
                                      , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$grth_rate_gam_max
                                      , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$grth_rate_gam_max
                                      , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$grth_rate_gam_max
                                      , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$grth_rate_gam_max
                                      , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$grth_rate_gam_max
                                      , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$grth_rate_gam_max
                                      , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$grth_rate_gam_max
                                      , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$grth_rate_gam_max
                                      , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$grth_rate_gam_max
                                      , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$grth_rate_gam_max
                                      , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$grth_rate_gam_max
                                      , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$grth_rate_gam_max
                                      , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$grth_rate_gam_max
                                      , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$grth_rate_gam_max
                                      , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$grth_rate_gam_max
)
tfps_lgr_mean_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                               , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$grth_rate_mean
                               , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$grth_rate_mean
                               , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$grth_rate_mean
                               , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$grth_rate_mean
                               , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$grth_rate_mean
                               , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$grth_rate_mean
                               , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$grth_rate_mean
                               , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$grth_rate_mean
                               , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$grth_rate_mean
                               , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$grth_rate_mean
                               , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$grth_rate_mean
                               , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$grth_rate_mean
                               , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$grth_rate_mean
                               , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$grth_rate_mean
                               , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$grth_rate_mean
                               , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$grth_rate_mean
                               , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$grth_rate_mean
                               , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$grth_rate_mean
                               , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$grth_rate_mean
                               , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$grth_rate_mean
                               , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$grth_rate_mean
                               , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$grth_rate_mean
                               , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$grth_rate_mean
                               , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$grth_rate_mean
)
tfps_lgr_simple_mean_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$grth_rate_simple_mean
                                , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$grth_rate_simple_mean
                                , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$grth_rate_simple_mean
                                , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$grth_rate_simple_mean
                                , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$grth_rate_simple_mean
                                , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$grth_rate_simple_mean
                                , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$grth_rate_simple_mean
                                , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$grth_rate_simple_mean
                                , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$grth_rate_simple_mean
                                , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$grth_rate_simple_mean
                                , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$grth_rate_simple_mean
                                , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$grth_rate_simple_mean
                                , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$grth_rate_simple_mean
                                , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$grth_rate_simple_mean
                                , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$grth_rate_simple_mean
                                , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$grth_rate_simple_mean
                                , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$grth_rate_simple_mean
                                , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$grth_rate_simple_mean
                                , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$grth_rate_simple_mean
                                , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$grth_rate_simple_mean
                                , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$grth_rate_simple_mean
                                , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$grth_rate_simple_mean
                                , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$grth_rate_simple_mean
                                , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$grth_rate_simple_mean
)
tfps_lgr_gam_mean_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$grth_rate_gam_mean
                                , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$grth_rate_gam_mean
                                , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$grth_rate_gam_mean
                                , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$grth_rate_gam_mean
                                , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$grth_rate_gam_mean
                                , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$grth_rate_gam_mean
                                , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$grth_rate_gam_mean
                                , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$grth_rate_gam_mean
                                , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$grth_rate_gam_mean
                                , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$grth_rate_gam_mean
                                , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$grth_rate_gam_mean
                                , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$grth_rate_gam_mean
                                , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$grth_rate_gam_mean
                                , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$grth_rate_gam_mean
                                , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$grth_rate_gam_mean
                                , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$grth_rate_gam_mean
                                , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$grth_rate_gam_mean
                                , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$grth_rate_gam_mean
                                , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$grth_rate_gam_mean
                                , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$grth_rate_gam_mean
                                , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$grth_rate_gam_mean
                                , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$grth_rate_gam_mean
                                , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$grth_rate_gam_mean
                                , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$grth_rate_gam_mean
)

tfps_lgr_simple_wtd_mean_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$grth_rate_simple_wtd_mean
                                , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$grth_rate_simple_wtd_mean
)

tfps_lgr_wtd_mean_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                    , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$grth_rate_wtd_mean
                                    , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$grth_rate_wtd_mean
                                    , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$grth_rate_wtd_mean
                                    , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$grth_rate_wtd_mean
                                    , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$grth_rate_wtd_mean
                                    , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$grth_rate_wtd_mean
                                    , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$grth_rate_wtd_mean
                                    , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$grth_rate_wtd_mean
                                    , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$grth_rate_wtd_mean
                                    , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$grth_rate_wtd_mean
                                    , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$grth_rate_wtd_mean
                                    , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$grth_rate_wtd_mean
                                    , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$grth_rate_wtd_mean
                                    , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$grth_rate_wtd_mean
                                    , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$grth_rate_wtd_mean
                                    , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$grth_rate_wtd_mean
                                    , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$grth_rate_wtd_mean
                                    , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$grth_rate_wtd_mean
                                    , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$grth_rate_wtd_mean
                                    , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$grth_rate_wtd_mean
                                    , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$grth_rate_wtd_mean
                                    , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$grth_rate_wtd_mean
                                    , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$grth_rate_wtd_mean
                                    , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$grth_rate_wtd_mean
)

tfps_lgr_gam_wtd_mean_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                    , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$grth_rate_gam_wtd_mean
                                    , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$grth_rate_gam_wtd_mean
)

#' Check all have the same number of dates
#length( tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date )
#length(  tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$growth_rate_var)
#length(  tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$growth_rate_var)
#length(  tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$growth_rate_var)
#length(  tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$growth_rate_var)
#length(  tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$growth_rate_var)
#length(  tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$growth_rate_var)
#length(  tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$growth_rate_var)
#length(  tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$growth_rate_var)
#length(  tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$growth_rate_var)
#length(  tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$growth_rate_var)
#length(  tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$growth_rate_var)
#length(  tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$growth_rate_var)
#length(  tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$growth_rate_var)
#length(  tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$growth_rate_var)
#length(  tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$growth_rate_var)
#length(  tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$growth_rate_var)
#length(  tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$growth_rate_var)
#length(  tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$growth_rate_var)
#' If there is a different number of rows
#a = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
#b = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$date
#a %ni% b
#tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date[214]

setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/analysis/2022_10_11 Leading indicator data - separate variables')
saveRDS( tfps_min_age_7_max_age_56_min_desc_20_growth_var_df , "mina07_maxa56_md020.rds" )
saveRDS( tfps_min_age_7_max_age_56_min_desc_50_growth_var_df , "mina07_maxa56_md050.rds" )
saveRDS( tfps_min_age_7_max_age_56_min_desc_100_growth_var_df , "mina07_maxa56_md100.rds" )
saveRDS( tfps_min_age_7_max_age_56_growth_var_df , "mina07_maxa56_mdperc.rds" )
saveRDS( tfps_min_age_7_max_age_84_min_desc_20_growth_var_df , "mina07_maxa84_md020.rds" )
saveRDS( tfps_min_age_7_max_age_84_min_desc_50_growth_var_df , "mina07_maxa84_md050.rds" )
saveRDS( tfps_min_age_7_max_age_84_min_desc_100_growth_var_df , "mina07_maxa84_md100.rds" )
saveRDS( tfps_min_age_7_max_age_84_growth_var_df , "mina07_maxa84_mdperc.rds" )

saveRDS( tfps_min_age_14_max_age_56_min_desc_20_growth_var_df , "mina14_maxa56_md020.rds" )
saveRDS( tfps_min_age_14_max_age_56_min_desc_50_growth_var_df , "mina14_maxa56_md050.rds" )
saveRDS( tfps_min_age_14_max_age_56_min_desc_100_growth_var_df , "mina14_maxa56_md100.rds" )
saveRDS( tfps_min_age_14_max_age_56_growth_var_df , "mina14_maxa56_mdperc.rds" )
saveRDS( tfps_min_age_14_max_age_84_min_desc_20_growth_var_df , "mina14_maxa84_md020.rds" )
saveRDS( tfps_min_age_14_max_age_84_min_desc_50_growth_var_df , "mina14_maxa84_md050.rds" )
saveRDS( tfps_min_age_14_max_age_84_min_desc_100_growth_var_df , "mina14_maxa84_md100.rds" )
saveRDS( tfps_min_age_14_max_age_84_growth_var_df , "mina14_maxa84_mdperc.rds" )

saveRDS( tfps_min_age_28_max_age_56_min_desc_20_growth_var_df , "mina28_maxa56_md020.rds" )
saveRDS( tfps_min_age_28_max_age_56_min_desc_50_growth_var_df , "mina28_maxa56_md050.rds" )
saveRDS( tfps_min_age_28_max_age_56_min_desc_100_growth_var_df , "mina28_maxa56_md100.rds" )
saveRDS( tfps_min_age_28_max_age_56_growth_var_df , "mina28_maxa56_mdperc.rds" )
saveRDS( tfps_min_age_28_max_age_84_min_desc_20_growth_var_df , "mina28_maxa84_md020.rds" )
saveRDS( tfps_min_age_28_max_age_84_min_desc_50_growth_var_df , "mina28_maxa84_md050.rds" )
saveRDS( tfps_min_age_28_max_age_84_min_desc_100_growth_var_df , "mina28_maxa84_md100.rds" )
saveRDS( tfps_min_age_28_max_age_84_growth_var_df , "mina28_maxa84_mdperc.rds" )


folder = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/analysis"
setwd( folder ) 
write.csv( tfps_vlgr_samp_df , file="tfps_vlgr_samp.csv")
write.csv( tfps_vlgr_simple_samp_df , file="tfps_vlgr_simple_samp.csv")
write.csv( tfps_vlgr_gam_samp_df , file="tfps_vlgr_gam_samp.csv")
write.csv( tfps_vlgr_pop_df , file="tfps_vlgr_pop.csv") 
write.csv( tfps_vlgr_simple_pop_df , file="tfps_vlgr_simple_pop.csv")
write.csv( tfps_vlgr_gam_pop_df , file="tfps_vlgr_gam_pop.csv")
write.csv( tfps_lgr_max_df , file="tfps_lgr_max.csv") 
write.csv( tfps_lgr_simple_max_df , file="tfps_lgr_simple_max.csv")
write.csv( tfps_lgr_gam_max_df , file="tfps_lgr_gam_max.csv")
write.csv( tfps_vlgr_mean_df , file="tfps_vlgr_mean.csv") 
write.csv( tfps_vlgr_simple_mean_df , file="tfps_vlgr_simple_mean.csv") 
write.csv( tfps_vlgr_gam_mean_df , file="tfps_vlgr_gam_mean.csv")
write.csv( tfps_vlgr_wtd_mean_df , file="tfps_vlgr_wtd_mean.csv")
write.csv( tfps_vlgr_simple_wtd_mean_df , file="tfps_vlgr_simple_wtd_mean.csv")
write.csv( tfps_vlgr_gam_wtd_mean_df , file="tfps_vlgr_gam_wtd_mean.csv")
write.csv( tfps_vlgr_wtd_df , file="tfps_vlgr_wtd.csv")
#write.csv( tfps_v_gam_lgr_df , file="tfps_v_gam_lgr.csv")
#write.csv( tfps_lead_ind_comp_df , file="tfps_lead_ind_comp.csv")
#write.csv( tfps_vlgr_mdperc_df , file="tfps_vlgr_mdperc.csv")


#'tfps_growth_rate_lists
tfps_growth_rate_lists = list(   "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_rate_list
                               , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_rate_list
                               , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_rate_list
                               , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_rate_list
                               , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_rate_list
                               , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_rate_list
                               , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_rate_list
                               , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_rate_list
                               , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_rate_list
                               , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_rate_list
                               , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_rate_list
                               , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_rate_list
                               , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_rate_list
                               , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_rate_list
                               , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_rate_list
                               , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_rate_list
                               , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_rate_list
                               , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_rate_list
                               , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_rate_list
                               , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_rate_list
                               , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_rate_list
                               , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_rate_list
                               , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_rate_list
                               , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_rate_list
                               )
tfps_growth_rate_gam_lists = list(   "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_rate_gam_list
                                 , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_rate_gam_list
                                 , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_rate_gam_list
                                 , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_rate_gam_list
                                 , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_rate_gam_list
                                 , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_rate_gam_list
                                 , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_rate_gam_list
                                 , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_rate_gam_list
                                 , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_rate_gam_list
                                 , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_rate_gam_list
                                 , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_rate_gam_list
                                 , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_rate_gam_list
                                 , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_rate_gam_list
                                 , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_rate_gam_list
                                 , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_rate_gam_list
                                 , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_rate_gam_list
                                 , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_rate_gam_list
                                 , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_rate_gam_list
                                 , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_rate_gam_list
                                 , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_rate_gam_list
                                 , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_rate_gam_list
                                 , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_rate_gam_list
                                 , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_rate_gam_list
                                 , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_rate_gam_list
) 
tfps_growth_rate_simple_lists = list(   "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_rate_simple_list
                                   , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_rate_simple_list
                                   , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_rate_simple_list
                                   , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_rate_simple_list
                                   , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_rate_simple_list
                                   , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_rate_simple_list
                                   , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_rate_simple_list
                                   , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_rate_simple_list
                                   , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_rate_simple_list
                                   , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_rate_simple_list
                                   , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_rate_simple_list
                                   , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_rate_simple_list
                                   , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_rate_simple_list
                                   , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_rate_simple_list
                                   , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_rate_simple_list
                                   , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_rate_simple_list
                                   , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_rate_simple_list
                                   , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_rate_simple_list
                                   , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_rate_simple_list
                                   , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_rate_simple_list
                                   , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_rate_simple_list
                                   , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_rate_simple_list
                                   , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_rate_simple_list
                                   , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_rate_simple_list
) 

tfps_clock_outlier_lists = list(   "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_clock_outlier_list
                                  , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_clock_outlier_list
                                  , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_clock_outlier_list
                                  , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_clock_outlier_list
                                  , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_clock_outlier_list
                                  , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_clock_outlier_list
                                  , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_clock_outlier_list
                                  , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_clock_outlier_list
                                  , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_clock_outlier_list
                                  , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_clock_outlier_list
                                  , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_clock_outlier_list
                                  , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_clock_outlier_list
                                  , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_clock_outlier_list
                                  , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_clock_outlier_list
                                  , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_clock_outlier_list
                                  , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_clock_outlier_list
                                  , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_clock_outlier_list
                                  , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_clock_outlier_list
                                  , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_clock_outlier_list
                                  , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_clock_outlier_list
                                  , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_clock_outlier_list
                                  , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_clock_outlier_list
                                  , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_clock_outlier_list
                                  , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_clock_outlier_list
)

tfps_cluster_size_lists = list(   "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_cluster_size_list
                                 , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_cluster_size_list
                                 , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_cluster_size_list
                                 , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_cluster_size_list
                                 , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_cluster_size_list
                                 , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_cluster_size_list
                                 , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_cluster_size_list
                                 , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_cluster_size_list
                                 , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_cluster_size_list
                                 , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_cluster_size_list
                                 , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_cluster_size_list
                                 , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_cluster_size_list
                                 , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_cluster_size_list
                                 , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_cluster_size_list
                                 , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_cluster_size_list
                                 , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_cluster_size_list
                                 , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_cluster_size_list
                                 , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_cluster_size_list
                                 , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_cluster_size_list
                                 , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_cluster_size_list
                                 , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_cluster_size_list
                                 , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_cluster_size_list
                                 , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_cluster_size_list
                                 , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_cluster_size_list
)

folder = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/analysis"
setwd( folder )
saveRDS( tfps_growth_rate_gam_lists , "tfps_growth_gam_lists.rds" )
saveRDS( tfps_growth_rate_simple_lists , "tfps_growth_simple_lists.rds" )
saveRDS( tfps_growth_rate_lists , "tfps_growth_lists.rds" )
saveRDS( tfps_cluster_size_lists , "tfps_cluster_size_lists.rds" )
saveRDS( tfps_clock_outlier_lists , "tfps_clock_outlier_lists.rds" )



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
library("grDevices")

#'MRC GIDA colour palette
patone_10 = rgb(138/256,121/256,103/256)  
patone_315 = rgb(33/256,103/256,126/256)  
patone_259 = rgb(106/256,59/256,119/256)
patone_228 = rgb(130/256,47/256,90/256)
patone_158 = rgb(208/256,114/256,50/256)
patone_382 = rgb(181/256,211/256,52/256)
patone_116 =rgb(255/256,219/256,0)

data_colour = c( patone_10 , patone_315 , patone_259 , patone_228 , patone_158 
                 , patone_382 , patone_116
                  ,"green", "black", "red", "purple", "blue", "orange"
                 , "darkgreen" , "yellow" , "aquamarine", "darkturquoise", "cadetblue" )

par(mfrow=c(1,1))
plot_df = tfps_vlgr_mdperc_df #tfps_vlgr_mdperc_df #tfps_vlgr_df[ , 1 ] # tfps_v_gam_lgr_df[ , 1 ] tfps_lgr_max_df[ , 1 ] tfps_vlgr_wtd_df[ , 1 ]
plot( plot_df[ , 1 ] 
      , plot_df[ , 2 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 4 ) #c( 0 , 4*18 ) 
      , col = data_colour
      , xaxp = c( min( plot_df[ , 1 ] ) , max( plot_df[ , 1 ] ) , 10 )
    )
abline(h=0,col="black")
for ( i in 3 : length( plot_df ) ){ 
  lines( plot_df[ , 1 ] 
         , plot_df[ , i ]
         , typ="l"
         , col=data_colour[ i-1 ])
}

points( hosp_df$date
        , hosp_df$cases/2000 +0.5
        #, typ="l"
        , col="grey")

#' rolling mean for vlgr / mean cluster size
plot_df = tfps_vlgr_max_df #tfps_vlgr_mdperc_df #tfps_vlgr_df[ , 1 ] # tfps_v_gam_lgr_df[ , 1 ] tfps_lgr_max_df[ , 1 ] tfps_vlgr_wtd_df[ , 1 ]
plot( plot_df[ 7:nrow(plot_df), 1 ] 
      , zoo::rollmean( plot_df[ , 2 ] , k = 7 )#, align = "right") 
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( -0 , 0.1 ) #c( 0 , 4*18 ) 
      , col = data_colour)
abline(h=0,col="black")
for ( i in 3 : length( plot_df ) ){ 
  lines( plot_df[ 7:nrow(plot_df), 1 ] 
         , zoo::rollmean( plot_df[ , i ] , k = 7 )
         , typ="l"
         , col=data_colour[ i-1 ])
}

points( hosp_df$date
        , hosp_df$cases/2000 +0.5
        #, typ="l"
        , col="grey")


#' Separate plots
plot( plot_df[ , 1 ]
      , plot_df[ , 2 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( -2 , 4*6 ) 
      , col = data_colour)

for ( i in 3 : length( plot_df ) ){
  message(i, " ",4*(i-2))
  lines( plot_df[ , 1 ]
         , plot_df[ , i ] + (4*(i-2))
         , typ="l"
         , col=data_colour[ i -1  ])
}
#lines(plot_df[ , 1 ], plot_df[ , 4 ]+70,typ="l")
points( hosp_df$date
        , hosp_df$cases/200
        , typ="l"
        , col="azure3")

#' min cluster age = 7, max cluseter age = 56, min descendants = 50
plot( tfps_vlgr_df[ , 1 ]
      , tfps_vlgr_df$mina07_maxa56_md050 # tfps_lgr_df[ , 4 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 4 ) 
      , col = data_colour[2])

points( hosp_df$date
        , hosp_df$cases/2000
        , typ="l"
        , col="azure3")
legend("topright"
       ,legend=c("min age = 7 days, max age = 56 days, min descendants = 50"
                           ,"hospitalisations (/2000)")
      ,lty=c(1,1)
      ,col=c(data_colour[2],"azure3")
      )

#' Other variance (growth rates) time series relative to min age=7, max age=56, min descendants=50
plot( tfps_vlgr_df[ , 1 ]
      , tfps_vlgr_df[ , 2 ] - tfps_vlgr_df[ , 3 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 4*18 ) 
      , col = data_colour)

for ( i in 3 : length( tfps_vlgr_df ) ){
  abline(h=0,col="grey")
  abline(h=4*(i-2),col="grey")
  lines( tfps_vlgr_df[ , 1 ]
         , tfps_vlgr_df[ , i ] - tfps_vlgr_df[ , 3 ] + (4*(i-2))
         , typ="l"
         , col=data_colour[ i -1 ])
}

points( hosp_df$date
        , hosp_df$cases/70
        , typ="l"
        , col="azure3")


#' Separate plots logistic growth rates vs GAM logistic growth rates
plot( tfps_vlgr_df[ , 1 ]
      , tfps_vlgr_df[ , 2 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 76)#4*18 ) 
      , col = "red" ) #data_colour)

for ( i in 3 : length( tfps_vlgr_df ) ){
  lines( tfps_vlgr_df[ , 1 ]
         , tfps_vlgr_df[ , i ] + (4*(i-2))
         , typ="l"
         , col= "red" ) # data_colour[ i ])
}

for ( i in 2 : length( tfps_v_gam_lgr_df ) ){
  lines( tfps_v_gam_lgr_df[ , 1 ]
         , tfps_v_gam_lgr_df[ , i ] + (4*(i-2))
         , typ="l"
         , col= "black" ) # data_colour[ i ])
}

points( hosp_df$date
        , hosp_df$cases/70
        , typ="l"
        , col="azure3" )

legend(  "topright"
       , legend=c(  "Logistic growth rate"
                  , "GAM logistic growth rate"
                  , "UK Covid-19 hospitalisations (1/70)")
                  ,lty=c(1,1,1)
        ,col = c( "red", "black", "azure3"))

#' Single example of lgr vs gam lgr . min cluster age = 7, max cluseter age = 56, min descendants = 100
plot( tfps_vlgr_df[ , 1 ]
      , tfps_vlgr_df[ , 2 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 4 ) 
      , col = "red")

lines( tfps_v_gam_lgr_df[ , 1 ]
      , tfps_v_gam_lgr_df[ , 2 ]
      , typ = "l" 
      , col = "black")

points( hosp_df$date
        , hosp_df$cases/2000
        , typ="l"
        , col="azure3")
legend("topright"
       ,legend=c(  "variance of logistic growth rate"
                 , "variance of GAM logistic growth rate"
                 , "hospitalisations (/2000)")
       ,lty=c(1,1)
       ,col=c("red","black","azure3")
)

#' Separate plots logistic growth rates vs logistic growth rates / mean cluster size
#' Need to normalise max values to 1 for each plot
plot( tfps_vlgr_df[ , 1 ]
      , tfps_vlgr_df[ , 2 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 76)#4*18 ) 
      , col = "red" ) #data_colour)

for ( i in 3 : length( tfps_vlgr_df ) ){
  lines( tfps_vlgr_df[ , 1 ]
         , tfps_vlgr_df[ , i ] + (4*(i-2))
         , typ="l"
         , col= "red" ) # data_colour[ i ])
}

for ( i in 2 : length( tfps_vlgr_df ) ){
  lines( tfps_vlgr_wtd_df[ , 1 ]
         , ( tfps_vlgr_wtd_df[ , i ] / max(tfps_vlgr_wtd_df[ , i ],na.rm=TRUE) * max(tfps_vlgr_df[ , i ],na.rm=TRUE) ) + (4*(i-2)) #' Normalised to max value of comparison time series + offset for plotting
         , typ="l"
         , col= "blue" ) # data_colour[ i ])
}

points( hosp_df$date
        , hosp_df$cases/70
        , typ="l"
        , col="azure3" )

legend(  "topright"
         , legend=c(  "Logistic growth rate"
                      , "Logistic growth rate / mean cluster size"
                      , "UK Covid-19 hospitalisations (1/70)")
         ,lty=c(1,1,1)
         ,col = c( "red", "blue", "azure3"))

#' Single example of lgr vs lgr/mean cluster size . min cluster age = 7, max cluseter age = 56, min descendants = 100
plot( tfps_vlgr_df[ , 1 ]
      , tfps_vlgr_df[ , 8 ]
      , typ = "l" 
      , xlab = "date"
      , ylab = "variance of cluster logistic growth rates"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 4 ) 
      , col = "red")

lines( tfps_vlgr_wtd_df[ , 1 ]
       , ( tfps_vlgr_wtd_df[ , 8 ] / max(tfps_vlgr_wtd_df[ , 8 ],na.rm=TRUE) * max(tfps_vlgr_df[ , 8 ],na.rm=TRUE) )  #' Normalised to max value of comparison time series + offset for plotting
       , typ="l"
       , col= "blue" ) # data_colour[ i ])

points( hosp_df$date
        , hosp_df$cases/2000
        , typ="l"
        , col="azure3")
legend("topright"
       ,legend=c(  "variance of logistic growth rate"
                   , "variance of logistic growth rate / mean cluster size"
                   , "hospitalisations (/2000)")
       ,lty=c(1,1)
       ,col=c("red","blue","azure3")
)

#' Leading indicator composite of multiple individual sets of variables. 
par(mfrow=c(1,1))

plot_x = tfps_vlgr_df$date
comp_1 = tfps_vlgr_df[ , 4 ]
comp_2 = tfps_vlgr_df[ , 5 ]
comp_3 = tfps_vlgr_df[ , 16 ]
comp_4 = tfps_lgr_max_df[ , 7 ]
comp_4[][is.infinite(comp_4[])] = NA
comp_a = (comp_1 + comp_2 + comp_3 + comp_4) / 4
comp_df = data.frame(    "comp1" = comp_1
                       , "comp2" = comp_2
                       , "comp3" = comp_3
                       , "comp4" = comp_4
                       )

comp_b = rowMeans( comp_df , na.rm = T )
comp_df = data.frame(  "date"  = plot_x
                     , "comp1" = comp_1
                     , "comp2" = comp_2
                     , "comp3" = comp_3
                     , "comp4" = comp_4
                     , "comp" = 
                     )

plot(plot_x,comp_a,typ="l")
lines(plot_x,comp_b,typ="l",col="red")

plot( plot_x
      , comp_b + 10
      , typ = "l" 
      , xlab = "date"
      , ylab = "leading indicator values"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 45 ) 
      , col = "black" #data_colour[1])
    )

lines( plot_x
       , comp_1 +15
       , typ="l"
       , col= "red" #data_colour[2] ) # data_colour[ i ])
      )

lines( plot_x
       , comp_2 + 20
       , typ="l"
       , col= "blue" #data_colour[3] ) # data_colour[ i ])
      )

lines( plot_x
       , comp_3 +25
       , typ="l"
       , col= "orange" #data_colour[5] ) # data_colour[ i ])
      )

lines( plot_x
       , comp_4 +30
       , typ="l"
       , col= "green" #data_colour[5] ) # data_colour[ i ])
        )
        
points( hosp_df$date
        , hosp_df$cases/500
        , typ="l"
        , col="grey")
legend("topright"
       ,legend=c(    "max logistic growth rate (min age = 7, max age = 84, min desc = 100)"
                   , "variance of logistic growth rate (min age = 28, max age = 56, min desc = 100)"
                   , "variance of logistic growth rate (min age = 7, max age = 84, min desc = 20)"
                   , "variance of logistic growth rate (min age = 7, max age = 56, min desc = 100)"
                   , "composite leading indicator"
                   , "hospitalisations (/500)")
       ,lty=c(1,1)
       #,col=c(data_colour[1],data_colour[2],data_colour[3],data_colour[5],"grey")
       ,col=c("green","orange","blue","red","black","grey")
        , cex = 0.75
)

tfps_lead_ind_comp_df = data.frame( "date" = plot_x
                                    , "composite" = comp_b
                                    )

#' Leading indicator composite (equal weight but individually normalised to 1) 
#' of multiple individual sets of variables. 
par(mfrow=c(1,1))
par(xaxt="s", yaxt="s")
par( mar = c( 5 , 5 , 4 , 4 ) ) #margins: bottom, left, top, right
plot_x = tfps_vlgr_df$date
comp_1 = tfps_vlgr_df[ , 4 ] / max(tfps_vlgr_df[ , 4 ],na.rm=T)
comp_2 = tfps_vlgr_df[ , 5 ] / max(tfps_vlgr_df[ , 5 ],na.rm=T)
comp_3 = tfps_vlgr_df[ , 16 ] / max(tfps_vlgr_df[ , 16 ],na.rm=T)
comp_4 = tfps_lgr_max_df[ , 7 ] / max( tfps_lgr_max_df[ , 7 ],na.rm=T)
comp_4[][is.infinite(comp_4[])] = NA
comp_a = (comp_1 + comp_2 + comp_3 + comp_4) / 4
comp_df = data.frame(    #"comp1" = comp_1
                          "comp2" = comp_2
                         , "comp3" = comp_3
                         , "comp4" = comp_4
)
comp_b = rowMeans( comp_df , na.rm = T )
comp_df = data.frame(  "date"  = plot_x
                       #, "comp1" = comp_1
                       , "comp2" = comp_2
                       , "comp3" = comp_3
                       , "comp4" = comp_4
                       , "comp" = comp_b
)

plot(plot_x,comp_a)#,typ="l")
lines(plot_x,comp_b,typ="l",col="red")

plot( plot_x
      , comp_b + 2
      , typ = "l" 
      , xlab = "date"
      , ylab = "leading indicator values"
      #, xlim = c( df_dates[ 150 , 1 ] , df_dates[ 288 , 1 ] )
      , ylim = c( 0 , 15 ) 
      , col = "black" #data_colour[1])
)

lines( plot_x
       , comp_1 +4
       , typ="l"
       , col= "red" #data_colour[2] ) # data_colour[ i ])
)

lines( plot_x
       , comp_2 + 6
       , typ="l"
       , col= "blue" #data_colour[3] ) # data_colour[ i ])
)

lines( plot_x
       , comp_3 +8
       , typ="l"
       , col= "orange" #data_colour[5] ) # data_colour[ i ])
)

lines( plot_x
       , comp_4 +10
       , typ="l"
       , col= "green" #data_colour[5] ) # data_colour[ i ])
)

points( hosp_df$date
        , hosp_df$cases/3000
        , typ="l"
        , col="grey")
legend("topright"
       ,legend=c(    "max logistic growth rate (min age = 7, max age = 84, min desc = 100)"
                     , "variance of logistic growth rate (min age = 28, max age = 56, min desc = 100)"
                     , "variance of logistic growth rate (min age = 7, max age = 84, min desc = 20)"
                     #, "variance of logistic growth rate (min age = 7, max age = 56, min desc = 100)"
                     , "composite leading indicator"
                     , "hospitalisations (/3000)")
       ,lty=c(1,1)
       #,col=c(data_colour[1],data_colour[2],data_colour[3],data_colour[5],"grey")
       ,col=c("green","orange","blue","black","grey")
       , cex = 0.75
)

tfps_lead_ind_comp_df = data.frame( "date" = plot_x
                                    , "composite" = comp_b
)

folder = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_09/Analysis"
setwd( folder )
write.csv( tfps_lead_ind_comp_df , file="tfps_lead_ind_comp.csv")

plot( hosp_df$date
      , hosp_df$cases
      , typ="l"
      , col="azure3",ylab = "hospitalisations",xlab="date")
points(plot_df[,1],replicate(nrow(plot_df),1),pch=3)
legend("topright", legend=c("hospitalisations","Tree dates"),col=c("azure3","black"),pch=c(NA,3),lty=c(1,NA))


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
