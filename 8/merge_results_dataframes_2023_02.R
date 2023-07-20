#' Takes data frames created on HPC using . These dataframes represent the 
#' generation of early warning signals (EWS) and assessment of True or False Positives
#' for different threshold levels. There are three types of data frame and a number of
#' files for each data frame relating to different values for the:
#'  - logistic growth rate threshold for replacing the sub-clusters with their parent cluster
#'  - the p-value threshold level for the cluster logistic growth rates
#'  - the type of leading indicator
#' 
#' This code takes all of the .rds files in a folder and combines them into three
#' dataframes depending on the names of the file

################### Create list of .rds files and check if any missing
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df/" )
file_list = list.files() ; file_list = subset( file_list , file_list != "desktop.ini" & file_list != "test" & file_list != "lgr_th_NA" )  
file_list_ews_results = as.data.frame( subset( file_list , grepl( "ews_results_df" , file_list ) ) )
file_list_wave_results = as.data.frame( subset( file_list , grepl( "wave_results_df" , file_list ) ) )
file_list_wave_results_analysis = as.data.frame( subset( file_list , grepl( "wave_results_analysis_df" , file_list ) ) )

#' Build list of expected file names
#min_age = c("07","14","28")
#max_age = c("56","84")
#min_desc = c("100","20","50","proportional")
p_val_th = c(0.01,0.05,10000)
parent_sub_lgr = c(0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1, 999)
leading_indicator = c(   "tfps_vlgr_samp"      , "tfps_vlgr_simple_samp"    , "tfps_vlgr_gam_samp" 
                       , "tfps_vlgr_pop"     , "tfps_vlgr_simple_pop"     , "tfps_vlgr_gam_pop"
                       , "tfps_vlgr_wtd"
                       , "tfps_lgr_max"      , "tfps_lgr_simple_max"      , "tfps_lgr_gam_max"
                       , "tfps_lgr_mean"     , "tfps_lgr_simple_mean"     , "tfps_lgr_gam_mean"
                       , "tfps_lgr_wtd_mean" , "tfps_lgr_simple_wtd_mean" , "tfps_lgr_gam_wtd_mean"
                       , "tfps_clock_outlier_max" , "tfps_clock_outlier_mean"
                       , "tfps_pango_lineage_max_lgr"
)

ews_results_check_list = data.frame()
wave_results_check_list = data.frame()
wave_results_analysis_check_list = data.frame()
for ( pval in p_val_th ){
  for ( psub in parent_sub_lgr ){
    for ( li in leading_indicator ){
      ews_results_new_row = as.data.frame( paste0(  "ews_results_df_lgr_th_" , psub 
                                        , "_pval_th_" , pval 
                                        , "_li_" , li , "_.rds"
                                      )
                              )
      ews_results_check_list = rbind( ews_results_check_list , ews_results_new_row )
      
      wave_results_new_row = as.data.frame( paste0(  "wave_results_df_lgr_th_" , psub 
                                                    , "_pval_th_" , pval 
                                                    , "_li_" , li , "_.rds"
      )
      )
      wave_results_check_list = rbind( wave_results_check_list , wave_results_new_row )
      
      wave_results_analysis_new_row = as.data.frame( paste0(  "wave_results_analysis_df_lgr_th_" , psub 
                                                     , "_pval_th_" , pval 
                                                     , "_li_" , li , "_.rds"
      )
      )
      wave_results_analysis_check_list = rbind( wave_results_analysis_check_list , wave_results_analysis_new_row )
    }
  }
}

library(sqldf)
missing_ews_results = sqldf('SELECT * FROM ews_results_check_list EXCEPT SELECT * FROM file_list_ews_results')
missing_wave_results = sqldf('SELECT * FROM wave_results_check_list EXCEPT SELECT * FROM file_list_wave_results')
missing_wave_results_analysis = sqldf('SELECT * FROM wave_results_analysis_check_list EXCEPT SELECT * FROM file_list_wave_results_analysis')

rm(  missing_ews_results    , missing_wave_results    , missing_wave_results_analysis 
   , ews_results_check_list , wave_results_check_list , wave_results_analysis_check_list
   , ews_results_new_row    , wave_results_new_row    , wave_results_analysis_new_row )

################# Bind files into single dataframes

setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df/" )
file_list = list.files() ; file_list = subset( file_list , file_list != "desktop.ini" & file_list != "test" & file_list != "lgr_th_NA" )  
file_list_ews_results = subset( file_list , grepl( "ews_results_df" , file_list ) )

#' Merge ews_results_df_XXX.rds
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df/" )
ews_results_df_list = lapply( file_list_ews_results , readRDS )
ews_results_df = data.table::rbindlist( ews_results_df_list )
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined/" )
saveRDS( ews_results_df , "ews_results_df.rds")
rm( file_list_ews_results , ews_results_df_list )
gc()

#' Merge wave_results_df_XXX.rds
#'** Uses a lot of memory - probably better to run this on the HPC**
file_list_wave_results = subset( file_list , grepl( "wave_results_df" , file_list ) )
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df/" )
wave_results_df_list = lapply( file_list_wave_results , readRDS ) # 20.9GB!!
wave_results_df = data.table::rbindlist( wave_results_df_list )
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined/" )
saveRDS( wave_results_df , "wave_results_df.rds")
rm( wave_results_df_list , file_list_wave_results )
gc()

#' Merge wave_results_analysis_results_df_XXX.rds
file_list_wave_results_analysis = subset( file_list , grepl( "wave_results_analysis_df" , file_list ) )
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df/" )
wave_results_analysis_df_list = lapply( file_list_wave_results_analysis , readRDS )
wave_results_analysis_df = data.table::rbindlist( wave_results_analysis_df_list )
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/EWS_df_combined/" )
saveRDS( wave_results_analysis_df , "wave_results_analysis_df.rds")
rm( wave_results_analysis_df_list , file_list_wave_results_analysis )
gc()
