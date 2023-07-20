#' Uses files produced on the high performance cluster using
#' "tfps_analysis_single_dir_filter_options_HPC_array.R"
#' and derivatives.
#' Input files are dataframes containing time series for various statistics 
#' compiled from TFP scans. Each dataframe is for a different set of variables 
#' (minimum cluster age, maximum cluster age, minimum number of descendants) 
#' e.g. 'tfps_min_age_7_max_age_56_min_desc_20_growth_var_df.rds'
#' 

#' Set up variables to cycle through to create folder names containing files produced using HPC
large_cluster_adjust_lgr_thresholds = c("060","065","070","075","080","085","090","095","100","false")
p_val_thresholds = c("no","001","005")

for ( lgr_th in large_cluster_adjust_lgr_thresholds ){
  for ( p_val_th in p_val_thresholds ){
    setwd( paste0("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_",lgr_th,"/p_val_filter_",p_val_th,"/dataframes_variables/") )
    file_list = list.files() ; file_list = subset( file_list , file_list != "desktop.ini" )
    
    for (i in 1 : length( file_list ) ) {
      temp_df <- readRDS( file_list[ i ] )  #temp_df <- readRDS( paste0(folder,file_list[ i ] ) )
      new_name = stringr::str_split( file_list[ i ] , pattern = ".rds")[[1]][1]
      assign( new_name , temp_df )
    }
    
    #' Compile statistics for different scan variables into single dataframes
    #' tfps_vlgr_samp_df (sample var instead of population var)
    #' tfps_vlgr_simple_samp_df (sample var instead of population var)
    #' tfps_vlgr_gam_samp_df (sample var instead of population var)
    #' tfps_vlgr_pop_df (population instead of sample variance)
    #' tfps_vlgr_simple_pop_df (sample var instead of population var)
    #' tfps_vlgr_gam_pop_df (sample var instead of population var)
    #' tfps_vlgr_wtd_df (sample var divided by mean cluster size)
    #' tfps_lgr_max_df
    #' tfps_lgr_simple_max_df
    #' tfps_lgr_gam_max_df
    #' tfps_lgr_mean_df
    #' tfps_lgr_simple_mean_df
    #' tfps_lgr_gam_mean_df
    #' tfps_lgr_wtd_mean_df
    #' tfps_lgr_simple_wtd_mean_df
    #' tfps_lgr_gam_wtd_mean_df
    #' tfps_clock_outlier_max_df
    #' tfps_clock_outlier_mean_df
    
    ##### tfps_vlgr_samp_df #####
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
    ##### tfps_vlgr_simple_samp_df #####
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
    ##### tfps_vlgr_gam_samp_df #####
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
    ##### tfps_vlgr_pop_df #####
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
    ##### tfps_vlgr_simple_pop_df #####
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
    ##### tfps_vlgr_gam_pop_df #####
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
    ##### tfps_vlgr_wtd_df #####
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
    ##### tfps_lgr_max_df #####
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
    ##### tfps_lgr_simple_max_df #####
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
    ##### tfps_lgr_gam_max_df #####
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
    ##### tfps_lgr_mean_df #####
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
    ##### tfps_lgr_simple_mean_df #####
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
    ##### tfps_lgr_gam_mean_df #####
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
    ##### tfps_lgr_simple_wtd_mean_df #####
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
    ##### tfps_lgr_wtd_mean_df #####
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
    ##### tfps_lgr_gam_wtd_mean_df #####
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
    ##### tfps_clock_outlier_max_df #####
    tfps_clock_outlier_max_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                             , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$clock_outlier_max
                                             , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$clock_outlier_max
                                             , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$clock_outlier_max
                                             , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$clock_outlier_max
                                             , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$clock_outlier_max
                                             , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$clock_outlier_max
                                             , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$clock_outlier_max
                                             , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$clock_outlier_max
                                             , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$clock_outlier_max
                                             , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$clock_outlier_max
                                             , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$clock_outlier_max
                                             , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$clock_outlier_max
                                             , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$clock_outlier_max
                                             , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$clock_outlier_max
                                             , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$clock_outlier_max
                                             , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$clock_outlier_max
                                             , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$clock_outlier_max
                                             , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$clock_outlier_max
                                             , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$clock_outlier_max
                                             , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$clock_outlier_max
                                             , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$clock_outlier_max
                                             , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$clock_outlier_max
                                             , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$clock_outlier_max
                                             , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$clock_outlier_max
    )
    ##### tfps_clock_outlier_mean_df #####
    tfps_clock_outlier_mean_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                              , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$clock_outlier_mean
                                              , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$clock_outlier_mean
                                              , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$clock_outlier_mean
                                              , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$clock_outlier_mean
                                              , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$clock_outlier_mean
                                              , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$clock_outlier_mean
                                              , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$clock_outlier_mean
                                              , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$clock_outlier_mean
                                              , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$clock_outlier_mean
                                              , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$clock_outlier_mean
                                              , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$clock_outlier_mean
                                              , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$clock_outlier_mean
                                              , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$clock_outlier_mean
                                              , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$clock_outlier_mean
                                              , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$clock_outlier_mean
                                              , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$clock_outlier_mean
                                              , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$clock_outlier_mean
                                              , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$clock_outlier_mean
                                              , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$clock_outlier_mean
                                              , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$clock_outlier_mean
                                              , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$clock_outlier_mean
                                              , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$clock_outlier_mean
                                              , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$clock_outlier_mean
                                              , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$clock_outlier_mean
    )
    ##### tfps_pango_lineage_max_df #####
    tfps_pango_lineage_max_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                             , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$pango_lineage_max
                                             , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$pango_lineage_max
                                             , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$pango_lineage_max
                                             , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$pango_lineage_max
                                             , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$pango_lineage_max
                                             , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$pango_lineage_max
                                             , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$pango_lineage_max
                                             , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$pango_lineage_max
                                             , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$pango_lineage_max
                                             , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$pango_lineage_max
                                             , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$pango_lineage_max
                                             , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$pango_lineage_max
                                             , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$pango_lineage_max
                                             , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$pango_lineage_max
                                             , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$pango_lineage_max
                                             , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$pango_lineage_max
                                             , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$pango_lineage_max
                                             , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$pango_lineage_max
                                             , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$pango_lineage_max
                                             , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$pango_lineage_max
                                             , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$pango_lineage_max
                                             , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$pango_lineage_max
                                             , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$pango_lineage_max
                                             , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$pango_lineage_max
    )
    ##### tfps_pango_lineage_max_lgr_df #####
    tfps_pango_lineage_max_lgr_df = data.frame(  "date" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
                                                 , "mina07_maxa56_md020" = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$pango_lineage_max_lgr
                                                 , "mina07_maxa56_md050" = tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$pango_lineage_max_lgr
                                                 , "mina07_maxa56_md100" = tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$pango_lineage_max_lgr
                                                 , "mina07_maxa84_md020" = tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$pango_lineage_max_lgr
                                                 , "mina07_maxa84_md050" = tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$pango_lineage_max_lgr
                                                 , "mina07_maxa84_md100" = tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$pango_lineage_max_lgr
                                                 , "mina14_maxa56_md020" = tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$pango_lineage_max_lgr
                                                 , "mina14_maxa56_md050" = tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$pango_lineage_max_lgr
                                                 , "mina14_maxa56_md100" = tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$pango_lineage_max_lgr
                                                 , "mina14_maxa84_md020" = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$pango_lineage_max_lgr
                                                 , "mina14_maxa84_md050" = tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$pango_lineage_max_lgr
                                                 , "mina14_maxa84_md100" = tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$pango_lineage_max_lgr
                                                 , "mina28_maxa56_md020" = tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$pango_lineage_max_lgr
                                                 , "mina28_maxa56_md050" = tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$pango_lineage_max_lgr
                                                 , "mina28_maxa56_md100" = tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$pango_lineage_max_lgr
                                                 , "mina28_maxa84_md020" = tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$pango_lineage_max_lgr
                                                 , "mina28_maxa84_md050" = tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$pango_lineage_max_lgr
                                                 , "mina28_maxa84_md100" = tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$pango_lineage_max_lgr
                                                 , "mina07_maxa56_mdperc" = tfps_min_age_7_max_age_56_growth_var_df$pango_lineage_max_lgr
                                                 , "mina07_maxa84_mdperc" = tfps_min_age_7_max_age_84_growth_var_df$pango_lineage_max_lgr
                                                 , "mina14_maxa56_mdperc" = tfps_min_age_14_max_age_56_growth_var_df$pango_lineage_max_lgr
                                                 , "mina14_maxa84_mdperc" = tfps_min_age_14_max_age_84_growth_var_df$pango_lineage_max_lgr
                                                 , "mina28_maxa56_mdperc" = tfps_min_age_28_max_age_56_growth_var_df$pango_lineage_max_lgr
                                                 , "mina28_maxa84_mdperc" = tfps_min_age_28_max_age_84_growth_var_df$pango_lineage_max_lgr
    )
    ##### Checks #####
    #' Check all have the same number of dates
    #length( tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date )
    #length(  tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date)
    #length(  tfps_min_age_7_max_age_56_min_desc_50_growth_var_df$date)
    #length(  tfps_min_age_7_max_age_56_min_desc_100_growth_var_df$date)
    #length(  tfps_min_age_7_max_age_84_min_desc_20_growth_var_df$date)
    #length(  tfps_min_age_7_max_age_84_min_desc_50_growth_var_df$date)
    #length(  tfps_min_age_7_max_age_84_min_desc_100_growth_var_df$date)
    #length(  tfps_min_age_14_max_age_56_min_desc_20_growth_var_df$date)
    #length(  tfps_min_age_14_max_age_56_min_desc_50_growth_var_df$date)
    #length(  tfps_min_age_14_max_age_56_min_desc_100_growth_var_df$date)
    #length(  tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$date)
    #length(  tfps_min_age_14_max_age_84_min_desc_50_growth_var_df$date)
    #length(  tfps_min_age_14_max_age_84_min_desc_100_growth_var_df$date)
    #length(  tfps_min_age_28_max_age_56_min_desc_20_growth_var_df$date)
    #length(  tfps_min_age_28_max_age_56_min_desc_50_growth_var_df$date)
    #length(  tfps_min_age_28_max_age_56_min_desc_100_growth_var_df$date)
    #length(  tfps_min_age_28_max_age_84_min_desc_20_growth_var_df$date)
    #length(  tfps_min_age_28_max_age_84_min_desc_50_growth_var_df$date)
    #length(  tfps_min_age_28_max_age_84_min_desc_100_growth_var_df$date)
    #' If there is a different number of rows
    #a = tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date
    #b = tfps_min_age_14_max_age_84_min_desc_20_growth_var_df$date
    #a %ni% b
    #tfps_min_age_7_max_age_56_min_desc_20_growth_var_df$date[214]
    
    #' These appear to be the same as the existing .rds files but saved with different names
    setwd( paste0( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_",lgr_th,"/p_val_filter_",p_val_th,"/" ) )
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
    
    setwd( paste0( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_",lgr_th,"/p_val_filter_",p_val_th,"/dataframes_statistics" ) )
    
    write.csv( tfps_vlgr_samp_df             , file="tfps_vlgr_samp.csv")
    write.csv( tfps_vlgr_simple_samp_df      , file="tfps_vlgr_simple_samp.csv")
    write.csv( tfps_vlgr_gam_samp_df         , file="tfps_vlgr_gam_samp.csv")
    write.csv( tfps_vlgr_pop_df              , file="tfps_vlgr_pop.csv") 
    write.csv( tfps_vlgr_simple_pop_df       , file="tfps_vlgr_simple_pop.csv")
    write.csv( tfps_vlgr_gam_pop_df          , file="tfps_vlgr_gam_pop.csv")
    write.csv( tfps_vlgr_wtd_df              , file="tfps_vlgr_wtd.csv") #' This is vlgr / mean cluster size - not sure this really makes sense as a leading indicator as var is higher at lower cluster size Chi-squared distribution
    write.csv( tfps_lgr_max_df               , file="tfps_lgr_max.csv") 
    write.csv( tfps_lgr_simple_max_df        , file="tfps_lgr_simple_max.csv")
    write.csv( tfps_lgr_gam_max_df           , file="tfps_lgr_gam_max.csv")
    write.csv( tfps_lgr_mean_df              , file="tfps_lgr_mean.csv")
    write.csv( tfps_lgr_simple_mean_df       , file="tfps_lgr_simple_mean.csv")
    write.csv( tfps_lgr_gam_mean_df          , file="tfps_lgr_gam_mean.csv")
    write.csv( tfps_lgr_wtd_mean_df          , file="tfps_lgr_wtd_mean.csv")
    write.csv( tfps_lgr_simple_wtd_mean_df   , file="tfps_lgr_simple_wtd_mean.csv")
    write.csv( tfps_lgr_gam_wtd_mean_df      , file="tfps_lgr_gam_wtd_mean.csv")
    write.csv( tfps_clock_outlier_max_df     , file="tfps_clock_outlier_max.csv")
    write.csv( tfps_clock_outlier_mean_df    , file="tfps_clock_outlier_mean.csv")
    write.csv( tfps_pango_lineage_max_df     , file="tfps_pango_lineage_max.csv")
    write.csv( tfps_pango_lineage_max_lgr_df , file="tfps_pango_lineage_max_lgr.csv")
    
    
    #write.csv( tfps_vlgr_mean_df , file="tfps_vlgr_mean.csv") 
    #write.csv( tfps_vlgr_simple_mean_df , file="tfps_vlgr_simple_mean.csv") 
    #write.csv( tfps_vlgr_gam_mean_df , file="tfps_vlgr_gam_mean.csv")
    #write.csv( tfps_vlgr_wtd_mean_df , file="tfps_vlgr_wtd_mean.csv")
    #write.csv( tfps_vlgr_simple_wtd_mean_df , file="tfps_vlgr_simple_wtd_mean.csv")
    #write.csv( tfps_vlgr_gam_wtd_mean_df , file="tfps_vlgr_gam_wtd_mean.csv")
    #write.csv( tfps_vlgr_wtd_df , file="tfps_vlgr_wtd.csv")
    #write.csv( tfps_v_gam_lgr_df , file="tfps_v_gam_lgr.csv")
    #write.csv( tfps_lead_ind_comp_df , file="tfps_lead_ind_comp.csv")
    #write.csv( tfps_vlgr_mdperc_df , file="tfps_vlgr_mdperc.csv")
    
    ###########################################
    
    #' Repeat process for the list

    setwd( paste0( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_",lgr_th,"/p_val_filter_",p_val_th,"/lists_variables" ) )

    file_list = list.files()
    
    for (i in 1 : length( file_list ) ) {
      new_name = stringr::str_split( file_list[ i ] , pattern = ".rds")[[1]][1]
      if ( new_name == "desktop.ini"){ next }
      temp_df = readRDS( file_list[ i ] )
      assign( new_name , temp_df )
    }
    
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
    
    #' Save list dataframes
    setwd( paste0( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_",lgr_th,"/p_val_filter_",p_val_th,"/lists_statistics" ) )
    
    saveRDS( tfps_growth_rate_gam_lists , "tfps_growth_gam_lists.rds" )
    saveRDS( tfps_growth_rate_simple_lists , "tfps_growth_simple_lists.rds" )
    saveRDS( tfps_growth_rate_lists , "tfps_growth_lists.rds" )
    saveRDS( tfps_cluster_size_lists , "tfps_cluster_size_lists.rds" )
    saveRDS( tfps_clock_outlier_lists , "tfps_clock_outlier_lists.rds" )
  }
}
