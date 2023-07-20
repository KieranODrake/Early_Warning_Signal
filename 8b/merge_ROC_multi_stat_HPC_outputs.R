#' Merges the data frames containing the ROC AUC values calculated for the TFPS leading indicators. 
#' There is a separate dataframe for each combination of:
#'  - logistic growth rate threshold for replacing the sub-clusters with their parent cluster
#'  - the p-value threshold level for the cluster logistic growth rates
#' 
#' This code takes all of the .rds files in a folder and combines them into three
#' dataframes depending on the names of the file

################### Create list of .rds files
#setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/ROC_multi_stat_HPC_5d/" )
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/ROC_multi_stat_HPC_10d/" )

file_list = list.files() ; file_list = subset( file_list , file_list != "desktop.ini" & file_list != "test" & file_list != "lgr_th_NA" )  
file_list_roc_multi_stat_rt = as.list( subset( file_list , grepl( "ROC_multi_stat_rt" , file_list ) ) )
file_list_roc_multi_stat_sd = as.list( subset( file_list , grepl( "ROC_multi_stat_sd" , file_list ) ) )

################# Bind files into single dataframes
#' Read RDS files
roc_multi_stat_rt_list = lapply( file_list_roc_multi_stat_rt , readRDS ) #' Rt critical transition as anchor for true positive period
roc_multi_stat_sd_list = lapply( file_list_roc_multi_stat_sd , readRDS ) #' Wave start(inflection) date as anchor for true positive period
#' Add columns for leading indicator variables
start.time <- Sys.time()
roc_multi_stat_rt_df = lapply( seq_along( roc_multi_stat_rt_list ) , 
                                     function( x , y, i ){
                                      new_col_mina = as.numeric( substr( x[[i]]$leading_indicator_type , 5 , 6 ) )
                                      new_col_maxa = as.numeric( substr( x[[i]]$leading_indicator_type , 12 , 13 ) )
                                      new_col_md = t( as.data.frame( lapply( seq_along( x[[i]]$leading_indicator_type ) , function( z , j ){
                                                                          if( substr( z[j] , 17 , 19 ) == "per" ){ 
                                                                            md = "perc"
                                                                          } else { 
                                                                            md = as.numeric( substr( z[j] , 17 , 19 ) ) 
                                                                          }
                                                                          return( md )  
                                                                        } , z = x[[i]]$leading_indicator_type )
                                                                ))
                                      start_char = t( as.data.frame( lapply( seq_along( x[[i]]$leading_indicator_type ) , function( z , j ){
                                                                              if( substr( z[j] , 17 , 19 ) == "per" ){ 
                                                                                start_char = 22
                                                                              } else { 
                                                                                start_char = 21
                                                                              }
                                                                              return( start_char )  
                                                                            } , z = x[[i]]$leading_indicator_type )
                                                     ))
                                      new_col_li = substr( x[[i]]$leading_indicator_type , start_char , nchar( x[[i]]$leading_indicator_type ) )
                      
                                      pslt = t( as.data.frame( lapply( seq_along( y ) , function( a , k ){
                                                              if( nchar( y[[k]] ) == max(nchar(y))-3 ){
                                                                pslt_lapply = substr( y[k] , 29 , 31 )
                                                              } else if ( nchar( y[[k]] ) == max(nchar(y))-2 ) { 
                                                                pslt_lapply = substr( y[k] , 29 , 31 )
                                                              } else if ( nchar( y[[k]] ) == max(nchar(y))-1 ) { 
                                                                pslt_lapply = substr( y[k] , 29 , 33 )
                                                              } else if ( nchar( y[[k]] ) == max(nchar(y)) ) { # max(nchar(y)) = 49 for 5d filenames and 53 for 10d filenames as '_10d' is appended
                                                                pslt_lapply = substr( y[k] , 29 , 33 )
                                                              }
                                                              return( pslt_lapply ) } , a = y )
                                                      ))
                                      pvt = t( as.data.frame( lapply( seq_along( y ) , function( a , k ){
                                                              if( nchar( y[[k]] ) == max(nchar(y))-3 ){
                                                                pvt_lapply  = substr( y[k] , 41 , 42 )
                                                              } else if ( nchar( y[[k]] ) == max(nchar(y))-2 ) { 
                                                                pvt_lapply  = substr( y[k] , 41 , 43 )
                                                              } else if ( nchar( y[[k]] ) == max(nchar(y))-1 ) { 
                                                                pvt_lapply  = substr( y[k] , 43 , 44 )
                                                              } else if ( nchar( y[[k]] ) == max(nchar(y)) ) { # max(nchar(y)) = 49 for 5d filenames and 53 for 10d filenames as '_10d' is appended
                                                                pvt_lapply  = substr( y[k] , 43 , 45 )
                                                              }
                                                              return( pvt_lapply ) } , a = y )
                                               ))
                        
                                      #' return dataframe with new columns added at the left and the original leading indicator labels removed
                                      return( cbind(   "leading_indicator" = as.data.frame( new_col_li )
                                                     , "min_age" = as.data.frame( new_col_mina )
                                                     , "max_age" = as.data.frame( new_col_maxa )
                                                     , "min_desc" = as.data.frame( new_col_md )
                                                     , "parent_sub_lgr_threshold" = as.data.frame( replicate( nrow( as.data.frame( new_col_li ) ) , pslt[i] ) )
                                                     , "pval_threshold" = as.data.frame( replicate( nrow( as.data.frame( new_col_li ) ) , pvt[i] ) )
                                                     , as.data.frame( x[[i]] )[-1] ) )
                                      }
                        , x = roc_multi_stat_rt_list , y = file_list_roc_multi_stat_rt )
#' Repeat for start date (sd)
roc_multi_stat_sd_df = lapply( seq_along( roc_multi_stat_sd_list ) , 
                               function( x , y, i ){
                                 new_col_mina = as.numeric( substr( x[[i]]$leading_indicator_type , 5 , 6 ) )
                                 new_col_maxa = as.numeric( substr( x[[i]]$leading_indicator_type , 12 , 13 ) )
                                 new_col_md = t( as.data.frame( lapply( seq_along( x[[i]]$leading_indicator_type ) , function( z , j ){
                                   if( substr( z[j] , 17 , 19 ) == "per" ){ 
                                     md = "perc"
                                   } else { 
                                     md = as.numeric( substr( z[j] , 17 , 19 ) ) 
                                   }
                                   return( md )  
                                 } , z = x[[i]]$leading_indicator_type )
                                 ))
                                 start_char = t( as.data.frame( lapply( seq_along( x[[i]]$leading_indicator_type ) , function( z , j ){
                                   if( substr( z[j] , 17 , 19 ) == "per" ){ 
                                     start_char = 22
                                   } else { 
                                     start_char = 21
                                   }
                                   return( start_char )  
                                 } , z = x[[i]]$leading_indicator_type )
                                 ))
                                 new_col_li = substr( x[[i]]$leading_indicator_type , start_char , nchar( x[[i]]$leading_indicator_type ) )
                                 
                                 pslt = t( as.data.frame( lapply( seq_along( y ) , function( a , k ){
                                   if( nchar( y[[k]] ) == max(nchar(y))-3 ){
                                     pslt_lapply = substr( y[k] , 29 , 31 )
                                   } else if ( nchar( y[[k]] ) == max(nchar(y))-2 ) { 
                                     pslt_lapply = substr( y[k] , 29 , 31 )
                                   } else if ( nchar( y[[k]] ) == max(nchar(y))-1 ) { 
                                     pslt_lapply = substr( y[k] , 29 , 33 )
                                   } else if ( nchar( y[[k]] ) == max(nchar(y)) ) { # max(nchar(y)) = 49 for 5d filenames and 53 for 10d filenames as '_10d' is appended
                                     pslt_lapply = substr( y[k] , 29 , 33 )
                                   }
                                   return( pslt_lapply ) } , a = y )
                                 ))
                                 pvt = t( as.data.frame( lapply( seq_along( y ) , function( a , k ){
                                   if( nchar( y[[k]] ) == max(nchar(y))-3 ){
                                     pvt_lapply  = substr( y[k] , 41 , 42 )
                                   } else if ( nchar( y[[k]] ) == max(nchar(y))-2 ) { 
                                     pvt_lapply  = substr( y[k] , 41 , 43 )
                                   } else if ( nchar( y[[k]] ) == max(nchar(y))-1 ) { 
                                     pvt_lapply  = substr( y[k] , 43 , 44 )
                                   } else if ( nchar( y[[k]] ) == max(nchar(y)) ) { # max(nchar(y)) = 49 for 5d filenames and 53 for 10d filenames as '_10d' is appended
                                     pvt_lapply  = substr( y[k] , 43 , 45 )
                                   }
                                   return( pvt_lapply ) } , a = y )
                                 ))
                                 
                                 #' return dataframe with new columns added at the left and the original leading indicator labels removed
                                 return( cbind(   "leading_indicator" = as.data.frame( new_col_li )
                                                  , "min_age" = as.data.frame( new_col_mina )
                                                  , "max_age" = as.data.frame( new_col_maxa )
                                                  , "min_desc" = as.data.frame( new_col_md )
                                                  , "parent_sub_lgr_threshold" = as.data.frame( replicate( nrow( as.data.frame( new_col_li ) ) , pslt[i] ) )
                                                  , "pval_threshold" = as.data.frame( replicate( nrow( as.data.frame( new_col_li ) ) , pvt[i] ) )
                                                  , as.data.frame( x[[i]] )[-1] ) )
                               }
                               , x = roc_multi_stat_sd_list , y = file_list_roc_multi_stat_sd )

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#' Bind dataframes
roc_multi_stat_rt_df = data.table::rbindlist( roc_multi_stat_rt_df ) 
roc_multi_stat_sd_df = data.table::rbindlist( roc_multi_stat_sd_df )

rm( roc_multi_stat_rt_list , roc_multi_stat_sd_list , file_list , file_list_roc_multi_stat_rt , file_list_roc_multi_stat_sd )
gc()

#' Change column names
colnames( roc_multi_stat_rt_df )[1:6] = c("leading_indicator","min_age","max_age","min_desc","parent_sub_lgr_threshold","pval_threshold")
colnames( roc_multi_stat_sd_df )[1:6] = c("leading_indicator","min_age","max_age","min_desc","parent_sub_lgr_threshold","pval_threshold")

#' Change format if necessary
unique( roc_multi_stat_rt_df$leading_indicator ) ; unique( roc_multi_stat_sd_df$leading_indicator )
unique( roc_multi_stat_rt_df$min_age ) ; unique( roc_multi_stat_sd_df$min_age )
unique( roc_multi_stat_rt_df$max_age ) ; unique( roc_multi_stat_sd_df$max_age )
unique( roc_multi_stat_rt_df$min_desc ) ; unique( roc_multi_stat_sd_df$min_desc )
unique( roc_multi_stat_rt_df$parent_sub_lgr_threshold ) ; unique( roc_multi_stat_sd_df$parent_sub_lgr_threshold )
unique( roc_multi_stat_rt_df$pval_threshold ) ; unique( roc_multi_stat_sd_df$pval_threshold )

#' Save dataframes
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/ROC_multi_stat_HPC_combined/" )
#saveRDS( roc_multi_stat_sd_df , "ROC_multi_stat_sd_tfps_combined_5d.rds")
#saveRDS( roc_multi_stat_rt_df , "ROC_multi_stat_rt_tfps_combined_5d.rds")
saveRDS( roc_multi_stat_sd_df , "ROC_multi_stat_sd_tfps_combined_10d.rds")
saveRDS( roc_multi_stat_rt_df , "ROC_multi_stat_rt_tfps_combined_10d.rds")

rm( AUC_sd_df , AUC_rt_df , file_list , file_list_auc_rt , file_list_auc_sd , auc_rt_list , auc_sd_list )
gc()

#' read in files
ROC_multi_stat_sd_tfps_combined_5d = readRDS( "ROC_multi_stat_sd_tfps_combined_5d.rds") ; ROC_multi_stat_sd_tfps_combined_5d = as.data.frame( ROC_multi_stat_sd_tfps_combined_5d )
ROC_multi_stat_rt_tfps_combined_5d = readRDS( "ROC_multi_stat_rt_tfps_combined_5d.rds") ; ROC_multi_stat_rt_tfps_combined_5d = as.data.frame( ROC_multi_stat_rt_tfps_combined_5d )

ROC_multi_stat_sd_tfps_combined_10d = readRDS( "ROC_multi_stat_sd_tfps_combined_10d.rds") ; ROC_multi_stat_sd_tfps_combined_10d = as.data.frame( ROC_multi_stat_sd_tfps_combined_10d )
ROC_multi_stat_rt_tfps_combined_10d = readRDS( "ROC_multi_stat_rt_tfps_combined_10d.rds") ; ROC_multi_stat_rt_tfps_combined_10d = as.data.frame( ROC_multi_stat_rt_tfps_combined_10d )

