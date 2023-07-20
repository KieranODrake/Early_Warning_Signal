#' Code to process/analyse output from tfp scanner. The product will be a data 
#' frame containing a time series of the variance of the cluster logistic growth
#' rates as well as other information. Several lists of lists will also be produced 
#' giving the individual values (such as cluster growth rates) for each scan date.
#' 
#' 1.1 - remove clusters without a recent sample/sequence
#' 1.2 - remove clusters with sub-clusters **removed on 13 Oct 22 as removing non-external clusters does what we actually want**
#' 1.3 - remove clusters that are not external
#' 1.4 - remove clusters with high p-values
#' 1.5 - remove clusters with overlapping tips
#' 1.6 - replace external clusters with parent clusters if parent growth rate is at least X% of max growth rate for sub_clusters
#' 2.1 - record logistic growth rates of remaining clusters
#' 2.2 - record cluster sizes of remaining clusters 
#' 2.3 - record clock outlier statistic
#' 3 - calculate variance of remaining cluster growth rates
#' 4 - output summary time series in data frames

library( magrittr ) 
library( lubridate ) 
#library( ggplot2 )
library( stringr )
library( zoo )

params_run <- commandArgs( trailingOnly = TRUE )

#' List of variables
min_age = c( 7 , 14 , 28 )
max_age = c( 56 , 84 )
min_desc = c( 20 , 50 , 100 )

#' Determine whether the minimum number of descendants (cluster size) is fixed 
#' or proportional to the number of samples in during the max-age period
min_desc_type = c( "fixed" , "proportional" )

#' Recreate file name suffix
counter = 0
fn_suffix = c()
for (i in 1 : length( min_age ) ){
  for (j in 1:length( max_age )){
    #' Next step depends on whether min descendants was fixed (and so consistently named in filenames)
    #' or is different depending on the date of the tree
    for ( mdt in min_desc_type ){
      if ( mdt == "fixed"){
        for (k in 1:length( min_desc )){
          counter = counter + 1
          fn_suffix[ counter ] = paste( "min_age_" , min_age[ i ] , "-max_age_" , max_age[ j ] , "-min_desc_" , min_desc[ k ], sep="" )
        }
      }
      if ( mdt == "proportional" ){
        counter = counter + 1
        fn_suffix[ counter ] = paste( "min_age_" , min_age[ i ] , "-max_age_" , max_age[ j ] , sep="" )
      }
    }
  }
}


rm( counter , i , j , k , mdt )

#' TFPScan file output of format 'scanner-2021-12-22-min_age_7-max_age_84-min_desc_20.rds'
#' 
#setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/outputs_perc/")
#setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/outputs/")
tfps_scan_folder = c( "mina07_maxa56_mind20" , "mina07_maxa56_mind50" , "mina07_maxa56_mind100" , "mina07_maxa56" , "mina07_maxa84_mind20" , "mina07_maxa84_mind50" , "mina07_maxa84_mind100" , "mina07_maxa84" , "mina14_maxa56_mind20" , "mina14_maxa56_mind50" , "mina14_maxa56_mind100" , "mina14_maxa56" , "mina14_maxa84_mind20" , "mina14_maxa84_mind50" , "mina14_maxa84_mind100" , "mina14_maxa84" , "mina28_maxa56_mind20" , "mina28_maxa56_mind50" , "mina28_maxa56_mind100" , "mina28_maxa56" , "mina28_maxa84_mind20" , "mina28_maxa84_mind50" , "mina28_maxa84_mind100" , "mina28_maxa84" )
n = as.numeric( params_run )
message( n )
message( params_run )
message( tfps_scan_folder )
message( tfps_scan_folder[ n ] )
message( paste( getwd() , "/" , tfps_scan_folder[ n ]  , sep = "" ))
setwd( paste( getwd() , "/" , tfps_scan_folder[ n ]  , sep = "" ) )
#fns = list.files( 'outputs', full.name=TRUE )
fns = list.files()

#' Select filters to apply to clusters (TRUE in each case means that filter will be applied)
#' 1.1 = remove clusters without a recent sample/sequence
extant = TRUE
#' 1.2 = remove clusters with sub-clusters
non_sub_clusters = FALSE
#' 1.3 = remove clusters that are not external
external = TRUE
#' 1.4 = replace sub-clusters with parent as long as the logistic growth rate is above threshold relative to sub-cluster maximum
#' Objective is to include more large clusters, which may produce better leading indicators for some waves
large_cluster_adjust = TRUE
parent_sub_lgr_threshold = 0.60 #' This should be varied
#' 1.5 = remove clusters with logistic growth p-values above a threshold
p_val_filter = TRUE
p_threshold = 0.01 #0.05 #10000
#' 1.6 = remove clusters with overlapping tips
non_overlap = TRUE

#' Cycle through sets of variables for analysis
#' Don't need for loop as running as array on High Performance Cluster
#for ( n in 1 : length( fn_suffix ) ) {
n = as.numeric( params_run )
  if ( grepl( "mind" , tfps_scan_folder[ n ] ) ) {
    #setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/outputs/")
    file_str_end = 51 #' For files where the number of minimum descendants is fixed
  } else {
    #setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_10/outputs_perc/")
    file_str_end = 39 #' For files where the number of minimum descendants is proportional to the number of sequences in the max cluster age period
  }
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
  grth_rate_var_pop_list <- c()        ; grth_rate_var_samp_list <- c()
  grth_rate_simple_var_pop_list <- c() ; grth_rate_simple_var_samp_list <- c()
  grth_rate_gam_var_pop_list <- c()    ; grth_rate_gam_var_samp_list <- c()
  #' Minimum for logistic growth rates, simple logistic growth rates and GAM logistic growth rates
  grth_rate_min_list <- c() ; grth_rate_simple_min_list <- c() ; grth_rate_gam_min_list <- c()
  #' Maximum for logistic growth rates, simple logistic growth rates and GAM logistic growth rates
  grth_rate_max_list <- c() ; grth_rate_simple_max_list <- c() ; grth_rate_gam_max_list <- c()
  #' Mean for logistic growth rates, simple logistic growth rates and GAM logistic growth rates
  grth_rate_mean_list <- c() ; grth_rate_simple_mean_list <- c() ; grth_rate_gam_mean_list <- c()
  #' Mean Weighted by cluster size for logistic growth rates, simple logistic growth rates and GAM logistic growth rates
  grth_rate_wtd_mean_list <- c() ; grth_rate_simple_wtd_mean_list <- c() ; grth_rate_gam_wtd_mean_list <- c()
  #' clock outlier stats - calculated on absolute values as sign is not important
  clock_outlier_max_list <- c( ) ; clock_outlier_mean_list <- c()
  #' For each pango lineage, find the maximum cluster (w/ contains most samples from given lineage) 
  # and return the maximum growth rate among all such clusters.
  pango_lineage_max_list <- c() ; pango_lineage_max_lgr_list <- c() 
  #' Number of clusters remaining after various filters
  n_clusters_raw_list <- c() ; n_clusters_na_grth_rate_list <- c() ; n_clusters_too_old_list <- c()
#  n_clusters_sub_list <- c()
  n_clusters_overlap_list <- c() ; n_clusters_remaining_list <- c() ; mean_cluster_size_list <- c()
  
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
    
      #' 1.1 - Remove clusters without a recent sample/sequence (14 days)
      if ( extant == TRUE ) {
        cut_off <-  max( tfps_output$most_recent_tip ) - 14
        
        tfps_output_filtered <- subset( tfps_output_filtered ,
                                        tfps_output_filtered$most_recent_tip >= cut_off
        )
        n_cluster_too_old <- n_cluster_raw - n_cluster_na_grth_rate - nrow( tfps_output_filtered )
      } else { n_cluster_too_old <- 0 }
      
      #'**do not use this part (2) as the objective is accomplished by section 4 filtering out internal (non-external) clusters**
      #' 1.2 - Remove clusters containing sub-clusters (and the sub-clusters)
      if ( non_sub_clusters == TRUE ){
        sub_clust <- intersect( tfps_output_filtered$cluster_id 
                                , tfps_output_filtered$parent_number )
        '%ni%' <- Negate("%in%") #' Define 'not in'
        tfps_output_filtered <- subset( tfps_output_filtered 
                                        , ( tfps_output_filtered$cluster_id %ni% sub_clust ) &
                                          ( tfps_output_filtered$parent_number %ni% sub_clust ) )
        n_cluster_sub <- n_cluster_raw - n_cluster_na_grth_rate - n_cluster_too_old - 
          nrow( tfps_output_filtered )
      } else { n_cluster_sub = 0 }
      
      #' 1.3 - Remove non-external clusters from dataset
      if ( external == TRUE ) {
        tfps_output_filtered = subset( tfps_output_filtered , external_cluster == TRUE )
      }
      
      #' 1.4 Replace some external clusters with the internal parent cluster where the internal parent growth rate is 
      #' at least X% of the maximum of the external sub-cluster growth rates
      #' Objective is to include more large clusters
      if ( large_cluster_adjust == TRUE & external == TRUE )  {
        #' Initialise dataframes
        rows_to_remove_list = data.frame()
        rows_to_add_list = data.frame()
        for ( parent in unique( tfps_output_filtered$parent_number ) ){
          #' look for external clusters with the same parent
          sub_cluster_df = subset( tfps_output_filtered , tfps_output_filtered$parent_number == parent )
          parent_cluster_df = subset( tfps_output , tfps_output$cluster_id == parent )
          if ( nrow( parent_cluster_df ) == 0 ){ next } #' Some parent clusters are not in tfps_output
          parent_lgr = parent_cluster_df$logistic_growth_rate
          
          #' if LGR of parent cluster is greater than x% of max( sub-cluster LGRs )... **NOTE THE EFFECT ON -VE LGR (smaller absolute value so less negative)** 
          if ( parent_lgr > ( parent_sub_lgr_threshold * max( sub_cluster_df$logistic_growth_rate ) ) ){
            #' ...then replace sub-clusters with parent cluster...
            #' ...but check the parent cluster against the 'extant' requirement and... 
            if ( ( extant == TRUE ) & ( parent_cluster_df$most_recent_tip >= cut_off ) ){
              #'...'LGR p-value' requirement before replacing
              if ( parent_cluster_df$logistic_growth_rate_p < p_threshold ){
                rows_to_remove_list = rbind( rows_to_remove_list , sub_cluster_df )
                rows_to_add_list = rbind( rows_to_add_list , parent_cluster_df )
              }
            }
          } else { next }#' if parent cluster doesn't meet criteria of LGR, extant and LGR p-value then assess next parent cluster
        }
        tfps_output_filtered = rbind( tfps_output_filtered , rows_to_add_list )
        tfps_output_filtered = tfps_output_filtered[ !tfps_output_filtered$cluster_id %in% rows_to_remove_list$cluster_id , ]
      }
      
      
      
      #' 1.5 - Remove clusters with high p-values for logistic growth rate above a threshold
      if ( p_val_filter == TRUE ){
        tfps_output_filtered = subset( tfps_output_filtered , logistic_growth_rate_p < p_threshold )
      }
      
      #' 1.6 - Remove clusters with overlapping tips
      if ( non_overlap == TRUE ){
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
      } else {
        tfps_output_filtered_overlap = tfps_output_filtered
        n_overlap_rows = 0
      }  
      
      #' 2.1 - Record set of logistic growth rates for clusters (after filtering) for each scan
      logistic_growth_rate_simple_list[[ i ]] <- tfps_output_filtered_overlap$simple_logistic_growth_rate
      logistic_growth_rate_gam_list[[ i ]] <- tfps_output_filtered_overlap$gam_logistic_growth_rate
      logistic_growth_rate_list[[ i ]] <- tfps_output_filtered_overlap$logistic_growth_rate
      
      #' 2.2 - Record set of cluster sizes for clusters (after filtering) for each scan
      cluster_size_list[[ i ]] <- tfps_output_filtered_overlap$cluster_size
      
      #' 2.3 - Record set of clock outlier statistics (after filtering) for each scan
      clock_outlier_list[[ i ]] <- tfps_output_filtered_overlap$clock_outlier
      
      #' 3 - Calculate variance of cluster logistic growth rates, weighted by cluster size
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
      
      if ( nrow( tfps_output_filtered_overlap ) != 0 ){
        #' Record clock outlier stats - calculated on absolute values as the sign is not important
        clock_outlier_max   <- max(  abs( tfps_output_filtered_overlap$clock_outlier ) , na.rm = TRUE)
        clock_outlier_mean  <- mean( abs( tfps_output_filtered_overlap$clock_outlier ) , na.rm = TRUE)
        
        #' Compute growth statistic
        #' For each pango lineage, find the maximum cluster (w/ contains most samples from given lineage) 
        #' and return the maximum growth rate among all such clusters. Repeat for all trees. 
        #' As per tfp_growth_stat_pango.R from Erik Volz on 7 Nov 2022
        #G1 <- sapply( ds, function(d){ #sapply( ds1, function(d){
        tfps_output_filtered_overlap$lineage <- sapply( strsplit( tfps_output_filtered_overlap$lineage, split = '\\|' ) , '[', 1) #d$lineage <- sapply( strsplit( d$lineage, split = '\\|' ) , '[', 1)
        lds <- split( tfps_output_filtered_overlap, tfps_output_filtered_overlap$lineage ) #lds <- split( d, d$lineage )
        lingr <- sapply( lds, function(d) d$simple_logistic_growth_rate[ which.max( d$cluster_size )] ) #lingr <- sapply( lds, function(d) d$simple_logistic_growth_rate[which.max(d$cluster_size)] )  
        k <- which.max( lingr )
        #print( length( lingr ) )
        #setNames( lingr[k], names( lds )[k] )
        pango_lineage_max <- names( lds )[k]
        pango_lineage_max_lgr <- data.frame(lingr[k])[1,1]
      }  else {
        clock_outlier_max     <- NA
        clock_outlier_mean    <- NA
        pango_lineage_max     <- NA
        pango_lineage_max_lgr <- NA
      }
      
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
      #n_clusters_sub_list            <- c( n_clusters_sub_list , n_cluster_sub )
      n_clusters_overlap_list        <- c( n_clusters_overlap_list , n_overlap_rows )
      n_clusters_remaining_list      <- c( n_clusters_remaining_list , number_clusters )
      mean_cluster_size_list         <- c( mean_cluster_size_list , mean( tfps_output_filtered_overlap$cluster_size ) )
      pango_lineage_max_list         <- c( pango_lineage_max_list , pango_lineage_max )
      pango_lineage_max_lgr_list     <- c( pango_lineage_max_lgr_list , pango_lineage_max_lgr )
      
      #end_time_single = Sys.time()
      #print( i )
      #print( end_time_single - start_time_single )
    }
    #end_time_total = Sys.time()
    #print( i )
    #print( end_time_total - start_time_total )
  
    #' 4 - output summary data for time series of scans/trees into data frame
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
                                  #, "n_clusters_sub"              <- n_clusters_sub_list
                                  , "n_clusters_overlap"          <- n_clusters_overlap_list
                                  , "n_clusters_remaining"        <- n_clusters_remaining_list
                                  , "mean_cluster_size"           <- mean_cluster_size_list
                                  , "pango_lineage_max"           <- pango_lineage_max_list
                                  , "pango_lineage_max_lgr"       <- pango_lineage_max_lgr_list 
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
                                     #, "n_clusters_sub"
                                     , "n_clusters_overlap"
                                     , "n_clusters_remaining"
                                     , "mean_cluster_size"
                                     , "pango_lineage_max"
                                     , "pango_lineage_max_lgr"
                                    )
    
    # Note that these dates are in numeric format and will need to be converted using as.Date( date_list[i], origin = "1970-01-01")
    names( logistic_growth_rate_list )        <- date_list 
    names( logistic_growth_rate_simple_list ) <- date_list 
    names( logistic_growth_rate_gam_list )    <- date_list 
    names( cluster_size_list )                <- date_list 
    names( clock_outlier_list )               <- date_list 
    
    #' Write to files: tfps_growth_var data frame and lists
    setwd("/rds/general/user/kdrake/home/tfps_2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_060/p_val_filter_001")

    var_values = gsub( "-" , "_" , fn_suffix[ n ] )
    saveRDS( tfps_growth_var                  , paste( "tfps_" , var_values , "_growth_var_df.rds"           , sep = "" ) )

    saveRDS( logistic_growth_rate_list        , paste( "tfps_" , var_values , "_growth_rate_list.rds"        , sep = "" ) )
    saveRDS( logistic_growth_rate_simple_list , paste( "tfps_" , var_values , "_growth_rate_simple_list.rds" , sep = "" ) )
    saveRDS( logistic_growth_rate_gam_list    , paste( "tfps_" , var_values , "_growth_rate_gam_list.rds"    , sep = "" ) )
    saveRDS( clock_outlier_list               , paste( "tfps_" , var_values , "_clock_outlier_list.rds"      , sep = "" ) )
    saveRDS( cluster_size_list                , paste( "tfps_" , var_values , "_cluster_size_list.rds"       , sep = "" ) )
