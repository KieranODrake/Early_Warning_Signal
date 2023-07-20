#' Compile dataframe containing information on highest growth variants and most frequent variants
#' for each set of scan variables (min age, max age, min descendants) and filter variable (LGR p-value) 
#' 
#' This is the version adapted for the Imperial High Performance Computing (HPC) environment using an array job

library( magrittr ) 
library( lubridate ) 
library( ggplot2 )
library( stringr )
library( zoo )
library( data.table )
'%!in%' = Negate('%in%')

#' Initialise data frame
top_variants_df = data.frame( "tree_date" = NA
                            , "TFPS_cluster_min_age" = NA
                            , "TFPS_cluster_max_age" = NA
                            , "TFPS_cluster_min_descendants" = NA
                            , "LGR_p_value_threshold" = NA
                            , "parent_sub_LGR_threshold" = NA
                            , "highest_freq_variant_1" = NA, "n_clusters_1" = NA
                            , "highest_freq_variant_2" = NA, "n_clusters_2" = NA
                            , "highest_freq_variant_3" = NA, "n_clusters_3" = NA
                            , "highest_freq_variant_4" = NA, "n_clusters_4" = NA
                            , "highest_freq_variant_5" = NA, "n_clusters_5" = NA
                            , "variant_1_LGR"=NA, "variant_1_LGR_freq"=NA,"variant_1_LGR_value"=NA,"variant_1_simple_LGR"=NA,"variant_1_simple_LGR_freq"=NA,"variant_1_simple_LGR_value"=NA,"variant_1_GAM_LGR"=NA,"variant_1_GAM_LGR_freq"=NA,"variant_1_GAM_LGR_value"=NA #' Information on variant with highest relevant growth rate
                            , "variant_2_LGR"=NA,"variant_2_LGR_freq"=NA,"variant_2_LGR_value"=NA,"variant_2_simple_LGR"=NA,"variant_2_simple_LGR_freq"=NA,"variant_2_simple_LGR_value"=NA,"variant_2_GAM_LGR"=NA,"variant_2_GAM_LGR_freq"=NA,"variant_2_GAM_LGR_value"=NA
                            , "variant_3_LGR"=NA,"variant_3_LGR_freq"=NA,"variant_3_LGR_value"=NA,"variant_3_simple_LGR"=NA,"variant_3_simple_LGR_freq"=NA,"variant_3_simple_LGR_value"=NA,"variant_3_GAM_LGR"=NA,"variant_3_GAM_LGR_freq"=NA,"variant_3_GAM_LGR_value"=NA
                            , "variant_4_LGR"=NA,"variant_4_LGR_freq"=NA,"variant_4_LGR_value"=NA,"variant_4_simple_LGR"=NA,"variant_4_simple_LGR_freq"=NA,"variant_4_simple_LGR_value"=NA,"variant_4_GAM_LGR"=NA,"variant_4_GAM_LGR_freq"=NA,"variant_4_GAM_LGR_value"=NA
                            , "variant_5_LGR"=NA,"variant_5_LGR_freq"=NA,"variant_5_LGR_value"=NA,"variant_5_simple_LGR"=NA,"variant_5_simple_LGR_freq"=NA,"variant_5_simple_LGR_value"=NA,"variant_5_GAM_LGR"=NA,"variant_5_GAM_LGR_freq"=NA,"variant_5_GAM_LGR_value"=NA
                            )

#############################################
#**HPC AMEND START**
params_run <- commandArgs( trailingOnly = TRUE )

#' List of variables used in TFP Scanner
#min_age = c( 7 , 14 , 28 )
#max_age = c( 56 , 84 )
#min_desc = c( 20 , 50 , 100 )
#min_desc = c( 20 , 50 , 100 , "proportional")
#**HPC AMEND END**

#' Select filters to apply to clusters (TRUE in each case means that filter will be applied)
#' These should be replicated from the version of tfps_analysis_single)dir_filter_options.R used to create the leading indicators being assessed for Early Warning Signals (EWS)
extant = TRUE #' remove clusters without a recent sample/sequence
non_sub_clusters = FALSE #' remove clusters with sub-clusters
external = TRUE #'remove clusters that are not external
large_cluster_adjust = TRUE #FALSE
#**HPC AMEND START**
lgr_threshold = c( 0.60 ) #0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 #* fix value in HPC
#  if ( large_cluster_adjust == TRUE ){
#    lgr_threshold = c( NA , seq( 0.60 , 1.00 , 0.05 ) ) #' Representing parent cluster LGR of between 60% and 100% of the maximum sub-cluster LGR
#  } else {
#    lgr_threshold = NA
#    lgr_th = NA
#  }

p_val_filter = TRUE #' remove clusters with logistic growth p-values above a threshold
#**HPC AMEND START**
p_threshold = c( 0.01 )#10000.00 , 0.05 , 0.01 fix value in HPC
#  if ( p_val_filter == TRUE ){
#    p_threshold = c( 10000.00 , 0.05 , 0.01 ) 
#  } else { p_threshold = c( 10000.00) }
#**HPC AMEND END**
non_overlap = TRUE #' remove clusters with overlapping tips

#' TFPS files are saved in folders named based on the scan variables, minimum cluster age (mina), maximum cluster age (maxa) and minimum descendants (mind)
#' The .pbs job file reads an array input file which will contain the index relating to tfp_scan_folder and determine the folder for reading in TFP scan files
tfps_scan_folder = c( "mina07_maxa56_mind20" , "mina07_maxa56_mind50" , "mina07_maxa56_mind100" , "mina07_maxa56" , "mina07_maxa84_mind20" , "mina07_maxa84_mind50" , "mina07_maxa84_mind100" , "mina07_maxa84" , "mina14_maxa56_mind20" , "mina14_maxa56_mind50" , "mina14_maxa56_mind100" , "mina14_maxa56" , "mina14_maxa84_mind20" , "mina14_maxa84_mind50" , "mina14_maxa84_mind100" , "mina14_maxa84" , "mina28_maxa56_mind20" , "mina28_maxa56_mind50" , "mina28_maxa56_mind100" , "mina28_maxa56" , "mina28_maxa84_mind20" , "mina28_maxa84_mind50" , "mina28_maxa84_mind100" , "mina28_maxa84" )
n = as.numeric( params_run )
message( n )
message( params_run )
message( tfps_scan_folder )
message( tfps_scan_folder[ n ] )
message( paste( getwd() , "/" , tfps_scan_folder[ n ]  , sep = "" ))
setwd( paste( getwd() , "/" , tfps_scan_folder[ n ]  , sep = "" ) )

#**HPC AMEND START**
file_list = list.files()
#' Cycle through sets of variables for analysis
#' Don't need for loop as running as array on High Performance Cluster
#for ( n in 1 : length( fn_suffix ) ) {
fn_suffix = tfps_scan_folder[ n ]
min_age = substr( tfps_scan_folder[ n ] , 5 , 6 )
max_age = substr( tfps_scan_folder[ n ] , 12 , 13 )
min_desc = ifelse( grepl( "mind" , tfps_scan_folder[ n ] ) ,  substr( tfps_scan_folder[ n ] , 19 , nchar(tfps_scan_folder[ n ]) ) , "proportional" )

#' Cycle through sets of scan variables and filter variables for analysis
  start_time_overall = Sys.time()
  #for ( a in 1 : length( min_age ) ){
  #  for ( b in 1 : length( max_age ) ){
  #    for ( c in 1 : length( min_desc ) ){
  #          #' Input files are outputs from the Transmission Fitness Polymorphism Scanner
  #          if ( min_desc[ c ] == "proportional" ){ #' when the minimum number of descendants in the cluster is set as proportional to the total number of sequences for that period
  #            setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_10/outputs_perc/")
  #            fns = list.files()
  #            fn_suffix = paste( "min_age_" , min_age[ a ] , "-max_age_" , max_age[ b ] , sep="" )
  #            file_str_end = 39
  #          } else { 
  #            setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_10/outputs/")
  #            fns = list.files()
  #            fn_suffix = paste( "min_age_" , min_age[ a ] , "-max_age_" , max_age[ b ] , "-min_desc_" , min_desc[ c ], sep="" )
  #            file_str_end = 51
  #          }
  #          #' This matches the middle part of the TFPScan filename 'scanner-2021-12-22-min_age_7-max_age_84-min_desc_20.rds'  
  #          file_list = subset( fns , substr( fns , 20 , file_str_end ) == fn_suffix )
  #          #' if file_list returns nothing then change the length of the character string to match
  #          if ( length( file_list ) == 0 ){
  #            file_list = subset( fns , substr( fns , 20 , file_str_end + 1 ) == fn_suffix )
  #          }
  #          if ( length( file_list ) == 0 ){
  #            file_list = subset( fns , substr( fns , 20 , file_str_end + 2 ) == fn_suffix )
  #          }
#**HPC AMEND end**
            #' Create list of scan dates
            scandates <- as.Date( regmatches( file_list , regexpr('(\\d{4}-\\d{2}-\\d{2})', file_list) ) )
      
            #' Loop through all TFP scans in list to obtain variants present in top 5 growth clusters
            #for ( f in 1: length( file_list ) ) {
              
              tfps_output = lapply( file_list , readRDS )
              names(tfps_output) <- scandates
              #' Remove NA LGRs
              tfps_output_filtered = lapply( tfps_output , function(x) { x[ !is.na( x$logistic_growth_rate ) , ] } )
              names(tfps_output_filtered) <- scandates
              
              for ( d in 1 : length( p_threshold ) ){
                #message( p_threshold[ d ] )
                for ( lgr_th in lgr_threshold ){
                  message( "Starting analysis of " , fn_suffix , ", p-val threshold: ", p_threshold[ d ], ", LGR threshold: " , lgr_th ," at ", Sys.time() )
                  
                start_time_individual = Sys.time()
                  
                  ###### Apply filters to clusters ############
                  #' 1.1 - Remove clusters without a recent sample/sequence (14 days)
                  if ( extant == TRUE ) {
                    tfps_output_filtered <- lapply( tfps_output_filtered 
                                                           , function( x ){ x[ x$most_recent_tip >= (max(x$most_recent_tip )-14) , ] } )
                  }
                  interval_time_1 = Sys.time()
                  message( "Filtered for extant clusters. Time taken: " , interval_time_1 - start_time_individual)
                
                  #' 1.3 - Remove non-external clusters from dataset
                  if ( external == TRUE ) {
                    tfps_output_filtered = lapply( tfps_output_filtered , function(x) { x[ x$external_cluster == TRUE , ] } )
                  }
                  interval_time_2 = Sys.time()
                  message( "Filtered for external clusters. Time taken: " , interval_time_2 - interval_time_1)
                  
                  #' 1.4 Replace some external clusters with the internal parent cluster where the internal parent growth rate is 
                  #' at least X% of the maximum of the external sub-cluster growth rates
                  #' Objective is to include more large clusters
                  #message("Analysing dataset " , fn_suffix , ", file " , f , ", p-val threshold: ", p_threshold[ d ], ", LGR threshold: " , lgr_th )
                  if ( is.na( lgr_th ) ){
                    #' Do nothing sd one of the parameters is to have no replacement of sub-clusters by parent clusters 
                  } else if ( large_cluster_adjust == TRUE & external == TRUE )  {
                    #' TEST https://stackoverflow.com/questions/9950144/access-lapply-index-names-inside-fun
                    #lapply( seq_along(tfps_output_filtered), function(y, z, n, m, i) {
                    #paste( y[[i]]$cluster_id[1] , z[[i]]$cluster_id[1] )
                    #}, y=tfps_output_filtered, z=tfps_output, n=names(tfps_output_filtered), m=names(tfps_output))
                    tfps_output_filtered <- lapply( seq_along( tfps_output_filtered ) , function( x , y, i ){
                                            cut_off <- max(x[[i]]$most_recent_tip )-14
                                            rows_to_remove_list = data.frame()
                                            rows_to_add_list = data.frame()
                                            for ( parent in unique( x[[i]]$parent_number ) ){
                                              #' look for external clusters with the same parent
                                              sub_cluster_df = subset( x[[i]] , x[[i]]$parent_number == parent )
                                              parent_cluster_df = subset( y[[i]] , y[[i]]$cluster_id == parent ); if ( nrow( parent_cluster_df ) == 0 ){ next } #' Some parent clusters are not in tfps_output
                                              parent_lgr = parent_cluster_df$logistic_growth_rate; if ( is.na( parent_lgr ) ){ next } #' Some parent clusters don't have a logistic growth rate
                                              
                                              #' if LGR of parent cluster is greater than x% of max( sub-cluster LGRs )... **NOTE THE EFFECT ON -VE LGR (smaller absolute value so less negative)** 
                                              if ( parent_lgr > ( lgr_th * max( sub_cluster_df$logistic_growth_rate ) ) ){
                                                #' ...then replace sub-clusters with parent cluster...
                                                #' ...but check the parent cluster against the 'extant' requirement and... 
                                                if ( ( extant == TRUE ) & ( parent_cluster_df$most_recent_tip >= cut_off ) ){
                                                  #'...'LGR p-value' requirement before replacing
                                                  if ( parent_cluster_df$logistic_growth_rate_p < p_threshold[ d ] ){
                                                    rows_to_remove_list = rbind( rows_to_remove_list , sub_cluster_df )
                                                    rows_to_add_list = rbind( rows_to_add_list , parent_cluster_df )
                                                  }
                                                }
                                              } else { next } #' if parent cluster doesn't meet criteria of LGR, extant and LGR p-value then assess next parent cluster
                                            }
                                            x[[i]] = rbind( x[[i]] , rows_to_add_list )
                                            #tfps_output_filtered = rbindlist( list( tfps_output_filtered , rows_to_add_list ) )
                                            x[[i]] = x[[i]][ !x[[i]]$cluster_id %in% rows_to_remove_list$cluster_id , ]
                                            #x[[i]] = x[[i]][ x[[i]]$cluster_id %!in% rows_to_remove_list$cluster_id , ]
                    } , x = tfps_output_filtered, y = tfps_output )  
                    names(tfps_output_filtered) <- scandates
                  }
                  interval_time_3 = Sys.time()
                  message( "Parent/sub-cluster replacement complete. Time taken: " , interval_time_3 - interval_time_2)
                  
                  #' 1.5 - Remove clusters with high p-values for logistic growth rate above a threshold
                  if ( p_val_filter == TRUE ){
                    tfps_output_filtered = lapply( tfps_output_filtered , function(x) { x[ x$logistic_growth_rate_p < p_threshold[ d ] , ] } )
                  }
                  interval_time_4 = Sys.time()
                  message( "Filtered for LGR p-value. Time taken: " , interval_time_4 - interval_time_3)
                  
                  #' 1.6 - Remove clusters with overlapping tips
                  if ( non_overlap == TRUE ){
                    
                    tfps_output_filtered_overlap = lapply( tfps_output_filtered , function(x) {
                      #' Check if there are any overlapping tips before trying to identify the overlapping clusters (which is more time consuming)
                      tips_all = data.frame( "tips" = c())
                      if ( nrow( x ) > 0 ){ #' sometimes all clusters are removed due to sub-clusters - so need to check
                        tips_all = data.frame( "tips" = c())
                        for ( j in 1 : nrow( x ) ) { 
                          tips_single = as.data.frame( strsplit( x$tips[ j ] , "|" , fixed = TRUE ) )
                          colnames( tips_single ) <- "tips"
                          tips_all <- as.data.frame( rbind( tips_all , tips_single ) )
                        }
                        if ( nrow( unique( tips_all ) ) < nrow( tips_all ) ){          overlap_assessment <- TRUE
                        } else if ( nrow( unique( tips_all ) ) == nrow( tips_all ) ) { overlap_assessment <- FALSE }
                      } else if ( nrow( x ) == 0 ){ overlap_assessment <- FALSE }
                      
                      #' Remove clusters that have overlapping tips (remove both clusters). 
                      #' Such clusters need to be removed as we only want to include 
                      #' non-overlapping clusters in our analysis.
                      if ( overlap_assessment == TRUE ){
                        overlap_rows <- c()
                        tip_freq <- data.frame( table( tips_all ) )
                        tips_overlap <- data.frame( subset( tip_freq , tip_freq$Freq > 1 )$tips )
                        colnames( tips_overlap ) <- "tips"
                        for ( k in 1 : nrow( x ) ) {
                          tips_k <- as.data.frame( strsplit( x$tips[ k ] , "|" , fixed = TRUE ) )
                          colnames( tips_k ) <- "tips"
                          n_overlaps = nrow( dplyr::intersect( tips_overlap , tips_k ) )
                          if ( n_overlaps > 0 ){
                            overlap_rows <- c( overlap_rows , k )
                          }
                        }
                        #n_overlap_rows <- length( overlap_rows )
                        tfps_output_filtered_overlap <- data.frame( x[ -overlap_rows , ] )
                      } else if ( overlap_assessment == FALSE ) { 
                        tfps_output_filtered_overlap = x
                        #n_overlap_rows = 0
                      }
                    })
                    
                  } else {
                    tfps_output_filtered_overlap = x
                    #n_overlap_rows = 0
                  }
                  
                  interval_time_5 = Sys.time()
                  message( "Overlapping clusters removed. Time taken: " , interval_time_5 - interval_time_4)
                  
                  #' Build dataframe of most frequent variants and highest growth variants
                  top_variants_list_final <- lapply( seq_along( tfps_output_filtered_overlap ) , function( x , m , i){ 
                            #' Add line of data to top_variants_df
                            #print(m)
                            #print(i)
                            #print(x[[i]])
                            new_row = top_variants_df[ 1 , ]
                            new_row$tree_date = m[i] #names( x[[i]] )[ i ] #tfps_date
                            new_row$TFPS_cluster_min_age = min_age #min_age[ a ]
                            new_row$TFPS_cluster_max_age = max_age #max_age[ b ]
                            new_row$TFPS_cluster_min_descendants = min_desc #min_desc[ c ] #'**could use a regexp to get the actual value**
                            new_row$LGR_p_value_threshold = p_threshold[ d ]
                            new_row$parent_sub_LGR_threshold = lgr_th
                            
                            x_df <- as.data.frame( x[[i]] )
                            names( x_df ) <- names( x[[1]] )
                            #print("made it")
                            
                            if( (nrow( x_df ) == 0) || (is.null( nrow( x_df ) )) ){ #' After filtering some scan outputs have no clusters remaining
                              top_variants_df = rbind ( top_variants_df , new_row )
                            } else {
                            
                              #' For the 5 largest growth rates, obtain the variant name and growth rate value
                              #' Logistic Growth Rate (LGR)
                              main_variants_df = data.frame()
                              main_variants_df_ordered = data.frame()
                              for (j in 1 : nrow( x_df ) ){
                                first_split = stringr::str_split( x_df[ j , ]$lineage_summary , "\n" )
                                main_variants_df[ j , 1 ] = stringr::str_split( first_split[[ 1 ]][ 3 ] , " " )[[ 1 ]][ 1 ]
                                main_variants_df[ j , 2 ] = stringr::str_split( first_split[[ 1 ]][ 3 ] , " " )[[ 1 ]][ 4 ]
                                main_variants_df[ j , 3 ] = x_df$logistic_growth_rate[ j ]
                                main_variants_df[ j , 4 ] = x_df$simple_logistic_growth_rate[ j ]
                                main_variants_df[ j , 5 ] = x_df$gam_logistic_growth_rate[ j ]
                                
                              }
                              names( main_variants_df ) = c("lineage_summary_highest_frequency_variant","variant_frequency","logistic_growth_rate","simple_logistic_growth_rate","GAM_logistic_growth_rate")
                              
                              #' Only record cluster variants where growth is > 0 as we are interested in the variants that are driving waves of infections
                              lgr_type = c("","simple_","GAM_")
                              
                              for ( item in lgr_type ){
                                lgr_col_name = paste( item , "logistic_growth_rate" , sep = "" )
                                main_variants_df_ordered <- subset( main_variants_df , !is.na( main_variants_df[ lgr_col_name ] ) )
                                main_variants_df_ordered <- data.table::setorderv( main_variants_df_ordered , cols = lgr_col_name , order = -1 )
                                
                                for ( n in 1 : 5 ){
                                  variable_column_name_1 = paste("variant_",n,"_",item,"LGR",sep="")
                                  variable_column_name_2 = paste("variant_",n,"_",item,"LGR_freq",sep="")
                                  variable_column_name_3 = paste("variant_",n,"_",item,"LGR_value",sep="")
                                  #' If is not NA & LGR is > 0
                                  if ( !is.na( main_variants_df_ordered[ n , lgr_col_name ] ) & ( main_variants_df_ordered[ n , lgr_col_name ] > 0 ) ){
                                      new_row[ , variable_column_name_1 ] = main_variants_df_ordered[ n , "lineage_summary_highest_frequency_variant" ]
                                      new_row[ , variable_column_name_2 ] = main_variants_df_ordered[ n , "variant_frequency" ]
                                      new_row[ , variable_column_name_3 ] = main_variants_df_ordered[ n , lgr_col_name ] 
                                  } else {
                                      new_row[ , variable_column_name_1 ] = ""
                                      new_row[ , variable_column_name_2 ] = ""
                                      new_row[ , variable_column_name_3 ] = ""
                                  }
                                }  
                              }
                              
                              #' Also add the variant that is the highest frequency in the largest number of clusters, but only with LGR > 0 
                              main_variants_dt = data.table::as.data.table( subset( main_variants_df , main_variants_df$logistic_growth_rate > 0 ) ) #' Same as line above but only use clusters with +ve LGR 
                              main_variants_dt = main_variants_dt[, .N ,by = lineage_summary_highest_frequency_variant ]
                              main_variants_dt = data.table::setorder( main_variants_dt , -N)
                              for ( n in 1 : 5 ){
                                high_freq_variant_col_name = paste( "highest_freq_variant_" , n , sep = "" )
                                new_row[ , high_freq_variant_col_name ] = main_variants_dt[ n , 1 ]
                                n_cluster_col_name = paste( "n_clusters_" , n , sep = "" )
                                new_row[ , n_cluster_col_name ]         = main_variants_dt[ n , 2 ]
                              }
                                
                              top_variants_df = rbindlist ( list( top_variants_df , new_row ) )
                            }
                            #' Remove the 1st row of NAs
                            #top_variants_df = lapply( top_variants_df , function(x){ x[-1,]})
                            
                  }, x = tfps_output_filtered_overlap , m = names(tfps_output_filtered_overlap) )
                  
                  interval_time_6 = Sys.time()
                  message( "Top variants dataframe built. Time taken: " , interval_time_6 - interval_time_5)
                  
                  top_variants_df_final = data.frame()
                  for ( i in 1:length(top_variants_list_final)){
                    temp = as.data.frame(top_variants_list_final[[i]])[-1,]
                    top_variants_df_final = rbindlist( list( top_variants_df_final , temp ) )
                  }
                  
                  end_time = Sys.time()
                  message( "Dataset analysed " , fn_suffix , ", p-val threshold: ", p_threshold[ d ], ", LGR threshold: " , lgr_th ,". Time taken - parameter: ", end_time - start_time_individual , ", - total so far:", end_time - start_time_overall)
                  #**HPC AMEND START**
                  #setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_12/analysis/top_variants/")
                  setwd("/rds/general/user/kdrake/home/tfps_2023_02/quantile_historical/analysis/top_variants/dataframes_output/")
                  #**HPC AMEND END**
                  saveRDS( top_variants_df_final , file = paste0("top_variants_",fn_suffix,"_p_val_",p_threshold[ d ],"_lgr_th_",lgr_th,".rds" ) )
                
                } #'lgr_threshold for loop end
                
            } #' p value for loop end
      #  }
     # }
    #}
  #}
#**HPC AMEND START**
#setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_12/analysis/top_variants/")
#saveRDS( top_variants_df , file = "top_variants.rds")
#write.csv( top_variants_df , file = "top_variants.csv")
#' When want to merge two top_variants_df dataframes created using different parameters
#tv_df = readRDS("top_variants_positive_large_df.rds")
#top_variants_combined = rbind( top_variants_df , tv_df )
#saveRDS( top_variants_combined , file = "top_variants_df.rds")
#write.csv( top_variants_combined , file = "top_variants_df.csv")
#**HPC AMEND END**