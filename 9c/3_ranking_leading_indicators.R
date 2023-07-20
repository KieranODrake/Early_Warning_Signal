#' Analysis of TFPS genomic and non-TFPS potential leading indicators for waves of SARS-CoV-2 in the UK

library(tibble)

#' Ranking the leading indicator parameter sets (genomic via TFPS and non-genomic) based on 
#' ROC AUC and normalised Matthews Correlation Coefficient (MCC)

#' Read in TFPS genomic leading indicator ROC stats
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/ROC_multi_stat_HPC_combined/" )
tfps_sd_5d = readRDS( "ROC_multi_stat_sd_tfps_combined_5d.rds") ; tfps_sd_5d = as.data.frame( tfps_sd_5d )
tfps_rt_5d = readRDS( "ROC_multi_stat_rt_tfps_combined_5d.rds") ; tfps_rt_5d = as.data.frame( tfps_rt_5d )
tfps_sd_10d = readRDS( "ROC_multi_stat_sd_tfps_combined_10d.rds") ; tfps_sd_10d = as.data.frame( tfps_sd_10d )
tfps_rt_10d = readRDS( "ROC_multi_stat_rt_tfps_combined_10d.rds") ; tfps_rt_10d = as.data.frame( tfps_rt_10d )

#' Read in non-genomic leading indicator stats
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/" )
non_tfps_sd_5d = readRDS("ROC_stats_sd_non_tfps_5d.rds")
non_tfps_rt_5d = readRDS("ROC_stats_rt_non_tfps_5d.rds")
non_tfps_sd_10d = readRDS("ROC_stats_sd_non_tfps_10d.rds")
non_tfps_rt_10d = readRDS("ROC_stats_rt_non_tfps_10d.rds")

#' Inspect results 
#'  Plots of all data - takes a long time to plot
par(mfrow=c(2,1))
plot(  x = tfps_rt_5d$w_all_ex_2_8_normMCC 
       , y = tfps_rt_5d$w_all_ex_2_8_ROC_AUC
       , col = "black"
       , xlim = c(0,1) , ylim = c(0,1)
       , xlab = "Normalised MCC (Alpha, Delta (1,2,3), BA.1)"
       , ylab = "ROC AUC (TPR vs FPR) (Alpha, Delta (1,2,3), BA.1)"
       , main = "'True' positive EWS window: -30<= t <= +5" )
points( x = tfps_sd_5d$w_all_ex_2_8_normMCC
        , y = tfps_sd_5d$w_all_ex_2_8_ROC_AUC
        , col="red")
abline( h = 0.7 , lty = 2 ) ; abline( v = 0.7 , lty = 2 )
legend(   "bottomright" 
          , legend = c( "Time window anchored to Rt critical transition" , "Time window anchored to wave start date" ) 
          , col = c("black","red") 
          , pch = c(16,16) )

plot(   x = tfps_rt_10d$w_all_ex_2_8_normMCC 
        , y = tfps_rt_10d$w_all_ex_2_8_ROC_AUC 
        , col="blue"
        , xlim = c(0,1) , ylim = c(0,1)
        , xlab = "Normalised MCC (Alpha, Delta (1,2,3), BA.1)"
        , ylab = "ROC AUC (TPR vs FPR) (Alpha, Delta (1,2,3), BA.1)"
        , main = "'True' positive EWS window: -30<= t <= +10" )
points( x = tfps_sd_10d$w_all_ex_2_8_normMCC 
        , y = tfps_sd_10d$w_all_ex_2_8_ROC_AUC 
        , col="green")
abline( h = 0.7 , lty = 2 ) ; abline( v = 0.7 , lty = 2 )
legend(   "bottomright" 
          , legend = c( "Time window anchored to Rt critical transition" , "Time window anchored to wave start date" ) 
          , col = c("blue","green") 
          , pch = c(16,16) )

#' Remove leading indicator sets with NaN for normMCC for waves 3,4,5,6,7 (i.e. those with sufficient data points)
tfps_sd_5d_filtered = subset( tfps_sd_5d , !is.nan( w3_normMCC ) & !is.nan( w4_normMCC ) & !is.nan( w5_normMCC ) & !is.nan( w6_normMCC ) & !is.nan( w7_normMCC ) )
tfps_rt_5d_filtered = subset( tfps_rt_5d , !is.nan( w3_normMCC ) & !is.nan( w4_normMCC ) & !is.nan( w5_normMCC ) & !is.nan( w6_normMCC ) & !is.nan( w7_normMCC ) )
tfps_sd_10d_filtered = subset( tfps_sd_10d , !is.nan( w3_normMCC ) & !is.nan( w4_normMCC ) & !is.nan( w5_normMCC ) & !is.nan( w6_normMCC ) & !is.nan( w7_normMCC ) )
tfps_rt_10d_filtered = subset( tfps_rt_10d , !is.nan( w3_normMCC ) & !is.nan( w4_normMCC ) & !is.nan( w5_normMCC ) & !is.nan( w6_normMCC ) & !is.nan( w7_normMCC ) )

non_tfps_sd_5d_filtered = subset( non_tfps_sd_5d , !is.nan( w3_normMCC ) & !is.nan( w4_normMCC ) & !is.nan( w5_normMCC ) & !is.nan( w6_normMCC ) & !is.nan( w7_normMCC ) )
non_tfps_rt_5d_filtered = subset( non_tfps_rt_5d , !is.nan( w3_normMCC ) & !is.nan( w4_normMCC ) & !is.nan( w5_normMCC ) & !is.nan( w6_normMCC ) & !is.nan( w7_normMCC ) )
non_tfps_sd_10d_filtered = subset( non_tfps_sd_10d , !is.nan( w3_normMCC ) & !is.nan( w4_normMCC ) & !is.nan( w5_normMCC ) & !is.nan( w6_normMCC ) & !is.nan( w7_normMCC ) )
non_tfps_rt_10d_filtered = subset( non_tfps_rt_10d , !is.nan( w3_normMCC ) & !is.nan( w4_normMCC ) & !is.nan( w5_normMCC ) & !is.nan( w6_normMCC ) & !is.nan( w7_normMCC ) )

print("Genomic, start date, +5 days:")      ; unique( tfps_sd_5d_filtered$leading_indicator )
print("Genomic, Rt, +5 days:")              ; unique( tfps_rt_5d_filtered$leading_indicator ) #' doesn't give any results when filter out NaNs for w3,w4,w5,w6,w7
print("Genomic, start date, +10 days:")     ; unique( tfps_sd_10d_filtered$leading_indicator )
print("Genomic, Rt, +10 days:")             ; unique( tfps_rt_10d_filtered$leading_indicator ) #' doesn't give any results when filter out NaNs for w3,w4,w5,w6,w7
print("Non-Genomic, start date, +5 days:")  ; unique( non_tfps_sd_5d_filtered$leading_indicator )
print("Non-Genomic, Rt, +5 days:")          ; unique( non_tfps_rt_5d_filtered$leading_indicator )
print("Non-Genomic, start date, +10 days:") ; unique( non_tfps_sd_10d_filtered$leading_indicator )
print("Non-Genomic, Rt, +10 days:")         ; unique( non_tfps_rt_10d_filtered$leading_indicator )

#' Add columns of analysis
#' Calculate minimum normMCC, arithmetic and geometric mean for waves 3 to 7
#'"Genomic, start date, +5 days:")
min_normMCC        = as.data.frame( apply( tfps_sd_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, min ) )
arith_mean_normMCC = as.data.frame( apply( tfps_sd_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, mean ) )
geom_mean_normMCC  = as.data.frame( apply( tfps_sd_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, function(x){ exp( mean( log(x) ) ) } ) )
tfps_sd_5d_filtered = cbind( tfps_sd_5d_filtered , min_normMCC , arith_mean_normMCC , geom_mean_normMCC )
names( tfps_sd_5d_filtered )[ c( ncol( tfps_sd_5d_filtered ) -2 , ncol( tfps_sd_5d_filtered ) -1 ,  ncol( tfps_sd_5d_filtered ) ) ] <- c( "min_normMCC_w3_to_w7" , "arith_mean_normMCC_w3_to_w7" , "geom_mean_normMCC_w3_to_w7" )
tfps_sd_5d_filtered = tfps_sd_5d_filtered[ with( tfps_sd_5d_filtered , order( w_all_ex_2_8_normMCC , decreasing = TRUE ) ) , ]
tfps_sd_5d_filtered_short = tfps_sd_5d_filtered[ , c( c(1:7) , 20,34,48,62,76,90,104, c( 146 : 150 ) ) ]
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/ROC_multi_stat_HPC_combined/filtered/" )
write.csv( tfps_sd_5d_filtered , "tfps_sd_5d_filtered.csv")
write.csv( tfps_sd_5d_filtered_short , "tfps_sd_5d_filtered_short.csv")
#'"Genomic, Rt critical transition, +5 days:")
min_normMCC        = as.data.frame( apply( tfps_rt_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, min ) )
arith_mean_normMCC = as.data.frame( apply( tfps_rt_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, mean ) )
geom_mean_normMCC  = as.data.frame( apply( tfps_rt_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, function(x){ exp( mean( log(x) ) ) } ) )
tfps_rt_5d_filtered = cbind( tfps_rt_5d_filtered , min_normMCC , arith_mean_normMCC , geom_mean_normMCC )
names( tfps_rt_5d_filtered )[ c( ncol( tfps_rt_5d_filtered ) -2 , ncol( tfps_rt_5d_filtered ) -1 ,  ncol( tfps_rt_5d_filtered ) ) ] <- c( "min_normMCC_w3_to_w7" , "arith_mean_normMCC_w3_to_w7" , "geom_mean_normMCC_w3_to_w7" )
tfps_rt_5d_filtered = tfps_rt_5d_filtered[ with( tfps_rt_5d_filtered , order( w_all_ex_2_8_normMCC , decreasing = TRUE ) ) , ]
tfps_rt_5d_filtered_short = tfps_rt_5d_filtered[ , c( c(1:7) , 20,34,48,62,76,90,104, c( 146 : 150 ) ) ]
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/ROC_multi_stat_HPC_combined/filtered/" )
write.csv( tfps_rt_5d_filtered , "tfps_rt_5d_filtered.csv")
write.csv( tfps_rt_5d_filtered_short , "tfps_rt_5d_filtered_short.csv")
#'"Genomic, start date, +10 days:")
min_normMCC        = as.data.frame( apply( tfps_sd_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, min ) )
arith_mean_normMCC = as.data.frame( apply( tfps_sd_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, mean ) )
geom_mean_normMCC  = as.data.frame( apply( tfps_sd_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, function(x){ exp( mean( log(x) ) ) } ) )
tfps_sd_10d_filtered = cbind( tfps_sd_10d_filtered , min_normMCC , arith_mean_normMCC , geom_mean_normMCC )
names( tfps_sd_10d_filtered )[ c( ncol( tfps_sd_10d_filtered ) -2 , ncol( tfps_sd_10d_filtered ) -1 ,  ncol( tfps_sd_10d_filtered ) ) ] <- c( "min_normMCC_w3_to_w7" , "arith_mean_normMCC_w3_to_w7" , "geom_mean_normMCC_w3_to_w7" )
tfps_sd_10d_filtered = tfps_sd_10d_filtered[ with( tfps_sd_10d_filtered , order( w_all_ex_2_8_normMCC , decreasing = TRUE ) ) , ]
tfps_sd_10d_filtered_short = tfps_sd_10d_filtered[ , c( c(1:7) , 20,34,48,62,76,90,104, c( 146 : 150 ) ) ]
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/ROC_multi_stat_HPC_combined/filtered/" )
write.csv( tfps_sd_10d_filtered , "tfps_sd_10d_filtered.csv")
write.csv( tfps_sd_10d_filtered_short , "tfps_sd_10d_filtered_short.csv")
#'"Genomic, Rt critical transition, +10 days:")
min_normMCC        = as.data.frame( apply( tfps_rt_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, min ) )
arith_mean_normMCC = as.data.frame( apply( tfps_rt_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, mean ) )
geom_mean_normMCC  = as.data.frame( apply( tfps_rt_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, function(x){ exp( mean( log(x) ) ) } ) )
tfps_rt_10d_filtered = cbind( tfps_rt_10d_filtered , min_normMCC , arith_mean_normMCC , geom_mean_normMCC )
names( tfps_rt_10d_filtered )[ c( ncol( tfps_rt_10d_filtered ) -2 , ncol( tfps_rt_10d_filtered ) -1 ,  ncol( tfps_rt_10d_filtered ) ) ] <- c( "min_normMCC_w3_to_w7" , "arith_mean_normMCC_w3_to_w7" , "geom_mean_normMCC_w3_to_w7" )
tfps_rt_10d_filtered = tfps_rt_10d_filtered[ with( tfps_rt_10d_filtered , order( w_all_ex_2_8_normMCC , decreasing = TRUE ) ) , ]
tfps_rt_10d_filtered_short = tfps_rt_10d_filtered[ , c( c(1:7) , 20,34,48,62,76,90,104, c( 146 : 150 ) ) ]
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/ROC_multi_stat_HPC_combined/filtered/" )
write.csv( tfps_rt_10d_filtered , "tfps_rt_10d_filtered.csv")
write.csv( tfps_rt_10d_filtered_short , "tfps_rt_10d_filtered_short.csv")
#'"Non-Genomic, start date, +5 days:")
min_normMCC        = as.data.frame( apply( non_tfps_sd_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, min ) )
arith_mean_normMCC = as.data.frame( apply( non_tfps_sd_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, mean ) )
geom_mean_normMCC  = as.data.frame( apply( non_tfps_sd_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, function(x){ exp( mean( log(x) ) ) } ) )
non_tfps_sd_5d_filtered = cbind( non_tfps_sd_5d_filtered , min_normMCC , arith_mean_normMCC , geom_mean_normMCC )
names( non_tfps_sd_5d_filtered )[ c( ncol( non_tfps_sd_5d_filtered ) -2 , ncol( non_tfps_sd_5d_filtered ) -1 ,  ncol( non_tfps_sd_5d_filtered ) ) ] <- c( "min_normMCC_w3_to_w7" , "arith_mean_normMCC_w3_to_w7" , "geom_mean_normMCC_w3_to_w7" )
non_tfps_sd_5d_filtered = non_tfps_sd_5d_filtered[ with( non_tfps_sd_5d_filtered , order( w_all_ex_2_8_normMCC , decreasing = TRUE ) ) , ]
non_tfps_sd_5d_filtered_short = non_tfps_sd_5d_filtered[ , c( c(1:2) , 15,29,43,57,71,85,99, c( 141 : 145 ) ) ]
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/filtered/" )
write.csv( non_tfps_sd_5d_filtered , "non_tfps_sd_5d_filtered.csv")
write.csv( non_tfps_sd_5d_filtered_short , "non_tfps_sd_5d_filtered_short.csv")
#'"Non-Genomic, Rt critical transition, +5 days:")
min_normMCC        = as.data.frame( apply( non_tfps_rt_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, min ) )
arith_mean_normMCC = as.data.frame( apply( non_tfps_rt_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, mean ) )
geom_mean_normMCC  = as.data.frame( apply( non_tfps_rt_5d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, function(x){ exp( mean( log(x) ) ) } ) )
non_tfps_rt_5d_filtered = cbind( non_tfps_rt_5d_filtered , min_normMCC , arith_mean_normMCC , geom_mean_normMCC )
names( non_tfps_rt_5d_filtered )[ c( ncol( non_tfps_rt_5d_filtered ) -2 , ncol( non_tfps_rt_5d_filtered ) -1 ,  ncol( non_tfps_rt_5d_filtered ) ) ] <- c( "min_normMCC_w3_to_w7" , "arith_mean_normMCC_w3_to_w7" , "geom_mean_normMCC_w3_to_w7" )
non_tfps_rt_5d_filtered = non_tfps_rt_5d_filtered[ with( non_tfps_rt_5d_filtered , order( w_all_ex_2_8_normMCC , decreasing = TRUE ) ) , ]
non_tfps_rt_5d_filtered_short = non_tfps_rt_5d_filtered[ , c( c(1:2) , 15,29,43,57,71,85,99, c( 141 : 145 ) ) ]
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/filtered/" )
write.csv( non_tfps_rt_5d_filtered , "non_tfps_rt_5d_filtered.csv")
write.csv( non_tfps_rt_5d_filtered_short , "non_tfps_rt_5d_filtered_short.csv")
#'"Non-Genomic, start date, +10 days:")
min_normMCC        = as.data.frame( apply( non_tfps_sd_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, min ) )
arith_mean_normMCC = as.data.frame( apply( non_tfps_sd_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, mean ) )
geom_mean_normMCC  = as.data.frame( apply( non_tfps_sd_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, function(x){ exp( mean( log(x) ) ) } ) )
non_tfps_sd_10d_filtered = cbind( non_tfps_sd_10d_filtered , min_normMCC , arith_mean_normMCC , geom_mean_normMCC )
names( non_tfps_sd_10d_filtered )[ c( ncol( non_tfps_sd_10d_filtered ) -2 , ncol( non_tfps_sd_10d_filtered ) -1 ,  ncol( non_tfps_sd_10d_filtered ) ) ] <- c( "min_normMCC_w3_to_w7" , "arith_mean_normMCC_w3_to_w7" , "geom_mean_normMCC_w3_to_w7" )
non_tfps_sd_10d_filtered = non_tfps_sd_10d_filtered[ with( non_tfps_sd_10d_filtered , order( w_all_ex_2_8_normMCC , decreasing = TRUE ) ) , ]
non_tfps_sd_10d_filtered_short = non_tfps_sd_10d_filtered[ , c( c(1:2) , 15,29,43,57,71,85,99, c( 141 : 145 ) ) ]
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/filtered/" )
write.csv( non_tfps_sd_10d_filtered , "non_tfps_sd_10d_filtered.csv")
write.csv( non_tfps_sd_10d_filtered_short , "non_tfps_sd_10d_filtered_short.csv")
#'"Non-Genomic, Rt critical transition, +10 days:")
min_normMCC        = as.data.frame( apply( non_tfps_rt_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, min ) )
arith_mean_normMCC = as.data.frame( apply( non_tfps_rt_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, mean ) )
geom_mean_normMCC  = as.data.frame( apply( non_tfps_rt_10d_filtered[ , c( "w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ], 1, function(x){ exp( mean( log(x) ) ) } ) )
non_tfps_rt_10d_filtered = cbind( non_tfps_rt_10d_filtered , min_normMCC , arith_mean_normMCC , geom_mean_normMCC )
names( non_tfps_rt_10d_filtered )[ c( ncol( non_tfps_rt_10d_filtered ) -2 , ncol( non_tfps_rt_10d_filtered ) -1 ,  ncol( non_tfps_rt_10d_filtered ) ) ] <- c( "min_normMCC_w3_to_w7" , "arith_mean_normMCC_w3_to_w7" , "geom_mean_normMCC_w3_to_w7" )
non_tfps_rt_10d_filtered = non_tfps_rt_10d_filtered[ with( non_tfps_rt_10d_filtered , order( w_all_ex_2_8_normMCC , decreasing = TRUE ) ) , ]
non_tfps_rt_10d_filtered_short = non_tfps_rt_10d_filtered[ , c( c(1:2) , 15,29,43,57,71,85,99, c( 141 : 145 ) ) ]
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/filtered/" )
write.csv( non_tfps_rt_10d_filtered , "non_tfps_rt_10d_filtered.csv")
write.csv( non_tfps_rt_10d_filtered_short , "non_tfps_rt_10d_filtered_short.csv")

rm( arith_mean_normMCC , geom_mean_normMCC , min_normMCC )

#' Merge tfps and non-tfps dataframes
#' First need to make sure they have the same columns (TFPS df have extra leading indicator variable columns) 
non_tfps_rt_10d_filtered = tibble::add_column( non_tfps_rt_10d_filtered , type = replicate(nrow(non_tfps_rt_10d_filtered),"Non-genomic")            , .before = "leading_indicator_type") #leading_indicator_type")
non_tfps_rt_10d_filtered = tibble::add_column( non_tfps_rt_10d_filtered , min_age = replicate(nrow(non_tfps_rt_10d_filtered),"NA")                  , .before = "ews_threshold") #leading_indicator_type")
non_tfps_rt_10d_filtered = tibble::add_column( non_tfps_rt_10d_filtered , max_age = replicate(nrow(non_tfps_rt_10d_filtered),"NA")                  , .before = "ews_threshold")
non_tfps_rt_10d_filtered = tibble::add_column( non_tfps_rt_10d_filtered , min_desc = replicate(nrow(non_tfps_rt_10d_filtered),"NA")                 , .before = "ews_threshold")
non_tfps_rt_10d_filtered = tibble::add_column( non_tfps_rt_10d_filtered , parent_sub_lgr_threshold = replicate(nrow(non_tfps_rt_10d_filtered),"NA") , .before = "ews_threshold")
non_tfps_rt_10d_filtered = tibble::add_column( non_tfps_rt_10d_filtered , pval_threshold = replicate(nrow(non_tfps_rt_10d_filtered),"NA")           , .before = "ews_threshold")
names(non_tfps_rt_10d_filtered)[2] <- "leading_indicator"

non_tfps_rt_5d_filtered = tibble::add_column( non_tfps_rt_5d_filtered , type = replicate(nrow(non_tfps_rt_5d_filtered),"Non-genomic")            , .before = "leading_indicator_type") #leading_indicator_type")
non_tfps_rt_5d_filtered = tibble::add_column( non_tfps_rt_5d_filtered , min_age = replicate(nrow(non_tfps_rt_5d_filtered),"NA")                  , .before = "ews_threshold") #leading_indicator_type")
non_tfps_rt_5d_filtered = tibble::add_column( non_tfps_rt_5d_filtered , max_age = replicate(nrow(non_tfps_rt_5d_filtered),"NA")                  , .before = "ews_threshold")
non_tfps_rt_5d_filtered = tibble::add_column( non_tfps_rt_5d_filtered , min_desc = replicate(nrow(non_tfps_rt_5d_filtered),"NA")                 , .before = "ews_threshold")
non_tfps_rt_5d_filtered = tibble::add_column( non_tfps_rt_5d_filtered , parent_sub_lgr_threshold = replicate(nrow(non_tfps_rt_5d_filtered),"NA") , .before = "ews_threshold")
non_tfps_rt_5d_filtered = tibble::add_column( non_tfps_rt_5d_filtered , pval_threshold = replicate(nrow(non_tfps_rt_5d_filtered),"NA")           , .before = "ews_threshold")
names(non_tfps_rt_5d_filtered)[2] <- "leading_indicator"

non_tfps_sd_10d_filtered = tibble::add_column( non_tfps_sd_10d_filtered , type = replicate(nrow(non_tfps_sd_10d_filtered),"Non-genomic")            , .before = "leading_indicator_type") #leading_indicator_type")
non_tfps_sd_10d_filtered = tibble::add_column( non_tfps_sd_10d_filtered , min_age = replicate(nrow(non_tfps_sd_10d_filtered),"NA")                  , .before = "ews_threshold") #leading_indicator_type")
non_tfps_sd_10d_filtered = tibble::add_column( non_tfps_sd_10d_filtered , max_age = replicate(nrow(non_tfps_sd_10d_filtered),"NA")                  , .before = "ews_threshold")
non_tfps_sd_10d_filtered = tibble::add_column( non_tfps_sd_10d_filtered , min_desc = replicate(nrow(non_tfps_sd_10d_filtered),"NA")                 , .before = "ews_threshold")
non_tfps_sd_10d_filtered = tibble::add_column( non_tfps_sd_10d_filtered , parent_sub_lgr_threshold = replicate(nrow(non_tfps_sd_10d_filtered),"NA") , .before = "ews_threshold")
non_tfps_sd_10d_filtered = tibble::add_column( non_tfps_sd_10d_filtered , pval_threshold = replicate(nrow(non_tfps_sd_10d_filtered),"NA")           , .before = "ews_threshold")
names(non_tfps_sd_10d_filtered)[2] <- "leading_indicator"

non_tfps_sd_5d_filtered = tibble::add_column( non_tfps_sd_5d_filtered , type = replicate(nrow(non_tfps_sd_5d_filtered),"Non-genomic")            , .before = "leading_indicator_type") #leading_indicator_type")
non_tfps_sd_5d_filtered = tibble::add_column( non_tfps_sd_5d_filtered , min_age = replicate(nrow(non_tfps_sd_5d_filtered),"NA")                  , .before = "ews_threshold") #leading_indicator_type")
non_tfps_sd_5d_filtered = tibble::add_column( non_tfps_sd_5d_filtered , max_age = replicate(nrow(non_tfps_sd_5d_filtered),"NA")                  , .before = "ews_threshold")
non_tfps_sd_5d_filtered = tibble::add_column( non_tfps_sd_5d_filtered , min_desc = replicate(nrow(non_tfps_sd_5d_filtered),"NA")                 , .before = "ews_threshold")
non_tfps_sd_5d_filtered = tibble::add_column( non_tfps_sd_5d_filtered , parent_sub_lgr_threshold = replicate(nrow(non_tfps_sd_5d_filtered),"NA") , .before = "ews_threshold")
non_tfps_sd_5d_filtered = tibble::add_column( non_tfps_sd_5d_filtered , pval_threshold = replicate(nrow(non_tfps_sd_5d_filtered),"NA")           , .before = "ews_threshold")
names(non_tfps_sd_5d_filtered)[2] <- "leading_indicator"

tfps_rt_10d_filtered = tibble::add_column( tfps_rt_10d_filtered , type = replicate(nrow(tfps_rt_10d_filtered),"Genomic") , .before = "leading_indicator") #leading_indicator_type")
tfps_rt_5d_filtered = tibble::add_column( tfps_rt_5d_filtered , type = replicate(nrow(tfps_rt_5d_filtered),"Genomic") , .before = "leading_indicator") #leading_indicator_type")
tfps_sd_10d_filtered = tibble::add_column( tfps_sd_10d_filtered , type = replicate(nrow(tfps_sd_10d_filtered),"Genomic") , .before = "leading_indicator") #leading_indicator_type")
tfps_sd_5d_filtered = tibble::add_column( tfps_sd_5d_filtered , type = replicate(nrow(tfps_sd_5d_filtered),"Genomic") , .before = "leading_indicator") #leading_indicator_type")

rt_10d_filtered = rbind( tfps_rt_10d_filtered , non_tfps_rt_10d_filtered )
rt_5d_filtered  = rbind( tfps_rt_5d_filtered  , non_tfps_rt_5d_filtered  )
sd_10d_filtered = rbind( tfps_sd_10d_filtered , non_tfps_sd_10d_filtered )
sd_5d_filtered  = rbind( tfps_sd_5d_filtered  , non_tfps_sd_5d_filtered  )

#' Explore relationship between different metrics
par(mfrow = c(3,4))
#par(mfrow = c(1,1))

tfps = tfps_sd_5d_filtered # tfps_sd_5d_filtered , tfps_sd_10d_filtered , tfps_rt_5d_filtered , tfps_rt_10d_filtered
non_tfps = non_tfps_sd_5d_filtered # non_tfps_sd_5d_filtered , non_tfps_sd_10d_filtered , non_tfps_rt_5d_filtered , non_tfps_rt_10d_filtered
tfps_list = list( tfps_sd_5d_filtered , tfps_sd_10d_filtered , tfps_rt_5d_filtered , tfps_rt_10d_filtered )
non_tfps_list = list( non_tfps_sd_5d_filtered , non_tfps_sd_10d_filtered , non_tfps_rt_5d_filtered , non_tfps_rt_10d_filtered )
#' Plot 4 x Min norm MCC vs Norm MCC
for (i in 1:4){
  tfps = as.data.frame( tfps_list[i] )
  non_tfps = as.data.frame( non_tfps_list[i] )
  plot( x = tfps$w_all_ex_2_8_normMCC
      , y = tfps$w_all_ex_2_8_ROC_AUC
      , xlim = c(0,1) , ylim = c(0,1)
      , xlab=""#, xlab = "Normalised MCC (Alpha, Delta (1,2,3), BA.1)" 
      , ylab = "ROC AUC" # ylab = "ROC AUC (Alpha, Delta (1,2,3), BA.1)"
      , cex.lab = 1.7
      , cex.axis = 1.7
      #, main = "Time period covering Alpha, Delta (1,2,3), Omicron BA.1 waves"
      ,cex.main = 1.7
      ,col="lightblue"
  )
  points( x = non_tfps$w_all_ex_2_8_normMCC 
      , y = non_tfps$w_all_ex_2_8_ROC_AUC
      , col = "black"  
      )
  abline(v=0.7,h=0.7,lty=2)
  #legend( "bottomright" 
  #      , legend=c("Genomic TFP Scanner","Non-genomic")
  #      , pch = c(16,16)
  #      , col = c("lightblue","black")
  #      , cex = 1.7)
}
#' Plot 4 x Min norm MCC vs Norm MCC
for (i in 1:4){
  tfps = as.data.frame( tfps_list[i] )
  non_tfps = as.data.frame( non_tfps_list[i] )
  plot( x = tfps$w_all_ex_2_8_normMCC 
    , y = tfps$min_normMCC_w3_to_w7 
    , xlim = c(0,1) , ylim = c(0,1)
    , xlab="" #, xlab = "Normalised MCC (Alpha, Delta (1,2,3), BA.1)" 
    , ylab = "Minimum normalised MCC" # ylab = "Minimum normalised MCC ((Alpha, Delta (1,2,3), BA.1)"
    , cex.lab = 1.7
    , cex.axis = 1.7
    ,col="green"
      )
  points( x = non_tfps$w_all_ex_2_8_normMCC 
        , y = non_tfps$min_normMCC_w3_to_w7
        , col = "black"  
  )
  abline(v=0.7,h=0.7,lty=2)
  #legend( "bottomright" 
  #      , legend=c("Genomic TFP Scanner","Non-genomic")
  #      , pch = c(16,16)
  #      , col = c("green","black")
  #      , cex = 1.7)
}

#' Plot 4 x Geometric mean of norm MCC vs Norm MCC
for (i in 1:4){
  tfps = as.data.frame( tfps_list[i] )
  non_tfps = as.data.frame( non_tfps_list[i] )
  plot( x = tfps$w_all_ex_2_8_normMCC
      , y = tfps$geom_mean_normMCC_w3_to_w7
      , xlim = c(0,1) , ylim = c(0,1)
      , xlab = "Normalised MCC (Alpha, Delta (1,2,3), BA.1)" 
      , ylab = "Geometric mean of normalised MCC" #, ylab = "Geometric mean of normalised MCC ((Alpha, Delta (1,2,3), BA.1)"
      , cex.lab = 1.7
      , cex.axis = 1.7
      ,col="yellow"
  )
  points( x = non_tfps$w_all_ex_2_8_normMCC 
        , y = non_tfps$geom_mean_normMCC_w3_to_w7
        , col = "black"  
  )
  abline(v=0.7,h=0.7,lty=2)
  #legend( "bottomright" 
  #      , legend=c("Genomic TFP Scanner","Non-genomic")
  #      , pch = c(16,16)
  #      , col = c("yellow","black")
  #      , cex = 1.7)
}

#' Violin plots by leading indicator https://plotly.com/r/violin/
library(plotly)

df <- rt_10d_filtered # read.csv("https://raw.githubusercontent.com/plotly/datasets/master/violin_data.csv")

x_var = ~leading_indicator #~type # ~leading_indicator
y_var = ~w_all_ex_2_8_normMCC# ~w_all_ex_2_normMCC #~w_all_normMCC #~w8_normMCC #~w_all_ex_2_8_normMCC #~w2_normMCC
y_title = "Normalised MCC" #"Normalised MCC (Alpha, Delta(1,2,3),BA.1)"

fig <- df %>%
  plot_ly(
    x = x_var, # ~leading_indicator, # ~day,
    y = y_var, # ~total_bill,
    split = x_var, #~leading_indicator, #~day,
    type = 'violin',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    )
  ) 

fig <- fig %>%
  layout(
    xaxis = list(
      title = "Leading indicator"#Day"
    ),
    yaxis = list(
      title = y_title, #"Total Bill",
      range = c(0,1),
      zeroline = F
    )
  )

fig

#fig <-fig + coord_flip()
#ggplotly(fig)

#' Violin plots for genomic and non-genomic leading indicators for start date / Rt critical transition and +5d / +10d
#' https://plotly.com/r/violin/
library(plotly)

rt_10d_filtered$time_window_anchor = replicate( nrow( rt_10d_filtered ) , "Rt+10")
rt_5d_filtered$time_window_anchor = replicate( nrow( rt_5d_filtered ) , "Rt+5")
sd_10d_filtered$time_window_anchor = replicate( nrow( sd_10d_filtered ) , "start_date+10")
sd_5d_filtered$time_window_anchor = replicate( nrow( sd_5d_filtered ) , "start_date+5")

combined_df <- rbind( rt_10d_filtered , rt_5d_filtered , sd_10d_filtered , sd_5d_filtered )
combined_df <- combined_df[ , c( 152, 1:151) ]
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/TFPS vs non-TFPS/" )
saveRDS( combined_df , "combined_df.rds")
write.csv( combined_df , "combined_df.csv")

#' Filter to reduce size and order by 
combined_filter_0.7_df = subset( combined_df , combined_df$w_all_ex_2_8_normMCC > 0.7 )
combined_filter_0.7_order_df = combined_filter_0.7_df[ with( combined_filter_0.7_df, order( w_all_ex_2_8_normMCC , decreasing = TRUE ) ) , ] 
View( combined_filter_0.7_order_df )
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/TFPS vs non-TFPS/" )
saveRDS( combined_filter_0.7_order_df , "combined_filter_07df.rds")
write.csv( combined_filter_0.7_order_df , "combined_filter_07_df.csv")

subset(combined_filter_0.7_order_df , combined_filter_0.7_order_df$time_window_anchor=="Rt+5" & combined_filter_0.7_order_df$type == "Non-genomic")
View( subset(combined_filter_0.7_order_df , combined_filter_0.7_order_df$time_window_anchor=="Rt+10" ) )

y_var = ~w2_normMCC# ~w_all_ex_2_normMCC #~w_all_normMCC #~w8_normMCC #~w_all_ex_2_8_normMCC #~w2_normMCC


fig <- combined_df %>%
  plot_ly(type = 'violin') 
fig <- fig %>%
  add_trace(
    x = ~time_window_anchor[combined_df$type == 'Genomic'], # ~day[df$smoker == 'Yes'],
    y = ~w_all_ex_2_8_normMCC[combined_df$type == 'Genomic'],
    legendgroup = 'Genomic',
    scalegroup = 'Genomic',
    name = 'Genomic',
    side = 'negative',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    ),
    color = I("blue")
  ) 
fig <- fig %>%
  add_trace(
    x = ~time_window_anchor[combined_df$type == 'Non-genomic'],
    y = ~w_all_ex_2_8_normMCC[combined_df$type == 'Non-genomic'],
    legendgroup = 'Non-genomic',
    scalegroup = 'Non-genomic',
    name = 'Non-genomic',
    side = 'positive',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    ),
    color = I("orange")
  ) 

fig <- fig %>%
  layout(
    xaxis = list(
      title = "",
      tickfont = list(size = 30)
    ),
    yaxis = list(
      title = "Normalised MCC",
      zeroline = F,
      range = c(0,1),
      tickfont = list(size = 30),
      titlefont = list(size = 30)
    ),
    legend = list(
      font = list(size=0),
      visible = F
    ),
    violingap = 0.3,
    violingroupgap = 0.3,
    violinmode = 'overlay'
  )

fig

#' Horizontal Violin plots by leading indicator https://plotly.com/r/violin/
library(plotly)

df <- rt_10d_filtered # read.csv("https://raw.githubusercontent.com/plotly/datasets/master/violin_data.csv")

y_var = ~leading_indicator #~type # ~leading_indicator
x_var = ~w_all_ex_2_8_normMCC# ~w_all_ex_2_normMCC #~w_all_normMCC #~w8_normMCC #~w_all_ex_2_8_normMCC #~w2_normMCC
y_title = "Leading indicator" #"Normalised MCC (Alpha, Delta(1,2,3),BA.1)"

fig <- df %>%
  plot_ly(
    x = x_var, # ~leading_indicator, # ~day,
    y = y_var, # ~total_bill,
    split = y_var, #~leading_indicator, #~day,
    type = 'violin',
    orientation="h",
    showlegend = FALSE,
    box = list( visible = T ),
    meanline = list( visible = T )
    #add_segments(x = 0.7)#, xend = 0.7)#, y = 0, yend = 10) #%>%
    #add_segments(x = 3, xend = 5, y = 5, yend = 5)
  ) 

vline <- function(x = 0, color = "black") {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color)
  )
}

fig <- fig %>%
  layout(
    xaxis = list(
      title = "Normalised MCC", #Day"
      range = c(0,1)
    ),
    yaxis = list(
      title = y_title, #"Total Bill",
      zeroline = F
    ),
    shapes = list( vline(0.7) )
  )

fig


# Filter for leading indicator sets with ROC AUC greater than 0.7 for period covering waves 3 to 7 
#'** actually probably better not to filter as will skew the percentile rankings**
data_filtered_AUC_07 = subset( data_filtered , w_all_ex_2_8_ROC_AUC > 0.7 )
points(data_filtered_AUC_07$w_all_ex_2_8_normMCC , data_filtered_AUC_07$w_all_ex_2_8_ROC_AUC , col="red")
#' Order by normalised MCC across period covering waves 3 to 7
data_filtered_AUC_07_order = data_filtered_AUC_07[ with( data_filtered_AUC_07 , order( w_all_ex_2_8_normMCC , decreasing = TRUE ) ) , ]
unique(data_filtered_AUC_07_order$leading_indicator)

#' Filter by norm MCC > 0.7 for period covering waves 3 to 7
#' Order 1st by normalised MCC across period covering waves 3 to 7
#'  secondly order by AUC, 
#'  3rdly by min Norm MCC and 
#'  4thly by geom mean norm MCC
combined_df_filter_normMCC_07 = subset( combined_df , w_all_ex_2_8_normMCC > 0.7 )
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/TFPS vs non-TFPS")
write.csv( combined_df_filter_normMCC_07 , "combined_df_filter_normMCC_07.csv")
combined_df_col_trim = combined_df[,c(1:9,22,36,50,64,78,92,106,117,121,148:152)]
write.csv( combined_df_col_trim , "combined_df_col_trim.csv")


#' Calculate mean, min across additional groups of waves (all with enough data: Alpha, Delta(1,2,3),BA.1 and genomic driven: Alpha, Delta(1), BA.1)
#' Add min and mean of genomic driven waves
min_normMCC_w3_4_7        = as.data.frame( apply( combined_df[ , c( "w3_normMCC","w4_normMCC","w7_normMCC") ], 1, min ) )
arith_mean_normMCC_w3_4_7 = as.data.frame( apply( combined_df[ , c( "w3_normMCC","w4_normMCC","w7_normMCC") ], 1, mean ) )
geom_mean_normMCC_w3_4_7  = as.data.frame( apply( combined_df[ , c( "w3_normMCC","w4_normMCC","w7_normMCC") ], 1, function(x){ exp( mean( log(x) ) ) } ) )
combined_df = cbind( combined_df , min_normMCC_w3_4_7 , arith_mean_normMCC_w3_4_7 , geom_mean_normMCC_w3_4_7 )
names( combined_df )[ c( ncol( combined_df ) -2 , ncol( combined_df ) -1 ,  ncol( combined_df ) ) ] <- c( "min_normMCC_w3_4_7" , "arith_mean_normMCC_w3_4_7" , "geom_mean_normMCC_w3_4_7" )
#' And add min and mean of non-genomic driven waves
min_normMCC_w5_6        = as.data.frame( apply( combined_df[ , c( "w5_normMCC","w6_normMCC") ], 1, min ) )
arith_mean_normMCC_w5_6 = as.data.frame( apply( combined_df[ , c( "w5_normMCC","w6_normMCC") ], 1, mean ) )
geom_mean_normMCC_w5_6  = as.data.frame( apply( combined_df[ , c( "w5_normMCC","w6_normMCC") ], 1, function(x){ exp( mean( log(x) ) ) } ) )
combined_df = cbind( combined_df , min_normMCC_w5_6 , arith_mean_normMCC_w5_6 , geom_mean_normMCC_w5_6 )
names( combined_df )[ c( ncol( combined_df ) -2 , ncol( combined_df ) -1 ,  ncol( combined_df ) ) ] <- c( "min_normMCC_w5_6" , "arith_mean_normMCC_w5_6" , "geom_mean_normMCC_w5_6" )

rm( min_normMCC_w3_4_7 , arith_mean_normMCC_w3_4_7 , geom_mean_normMCC_w3_4_7 , min_normMCC_w5_6 , arith_mean_normMCC_w5_6 , geom_mean_normMCC_w5_6  )

setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/filtered/" )
write.csv( non_tfps_rt_10d_filtered , "non_tfps_rt_10d_filtered.csv")
write.csv( non_tfps_rt_10d_filtered_short , "non_tfps_rt_10d_filtered_short.csv")

#' Calculate percentile ranking
#' All waves with sufficient data: Alpha, Delta(1,2,3),BA.1
#combined_df$w3_7_normMCC_perc_rank = ( rank( combined_df$w_all_ex_2_8_normMCC ) / length( combined_df$w_all_ex_2_8_normMCC ) )* 100
combined_df$w3_4_5_6_7_min_normMCC_perc_rank = ( rank( combined_df$min_normMCC_w3_to_w7 ) / length( combined_df$min_normMCC_w3_to_w7 ) )* 100
combined_df$w3_4_5_6_7_geom_mean_normMCC_perc_rank = ( rank( combined_df$geom_mean_normMCC_w3_to_w7 ) / length( combined_df$geom_mean_normMCC_w3_to_w7 ) )* 100
#combined_df$w3_7_ROC_AUC_perc_rank = ( rank( combined_df$w_all_ex_2_8_ROC_AUC ) / length( combined_df$w_all_ex_2_8_ROC_AUC ) )* 100
combined_df$w3_4_5_6_7_perc_rank_mean = rowMeans(  data.frame( combined_df$w3_4_5_6_7_min_normMCC_perc_rank
                                                             , combined_df$w3_4_5_6_7_geom_mean_normMCC_perc_rank
                                                             )
                                                )
#' Genomic-driven waves with sufficient data: Alpha, Delta(1),BA.1
combined_df$w3_4_7_min_normMCC_perc_rank = ( rank( combined_df$min_normMCC_w3_4_7 ) / length( combined_df$min_normMCC_w3_4_7 ) )* 100
combined_df$w3_4_7_geom_mean_normMCC_perc_rank = ( rank( combined_df$geom_mean_normMCC_w3_4_7 ) / length( combined_df$geom_mean_normMCC_w3_4_7 ) )* 100
combined_df$w3_4_7_perc_rank_mean = rowMeans(  data.frame( combined_df$w3_4_7_min_normMCC_perc_rank
                                                               , combined_df$w3_4_7_geom_mean_normMCC_perc_rank
                                                          )
                                            )
#' Non-Genomic-driven waves with sufficient data: Delta(2,3)
combined_df$w5_6_min_normMCC_perc_rank = ( rank( combined_df$min_normMCC_w5_6 ) / length( combined_df$min_normMCC_w5_6 ) )* 100
combined_df$w5_6_geom_mean_normMCC_perc_rank = ( rank( combined_df$geom_mean_normMCC_w5_6 ) / length( combined_df$geom_mean_normMCC_w5_6 ) )* 100
combined_df$w5_6_perc_rank_mean = rowMeans(  data.frame( combined_df$w5_6_min_normMCC_perc_rank
                                                           , combined_df$w5_6_geom_mean_normMCC_perc_rank
                                                        )
                                          )
#' save data frame
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/TFPS vs non-TFPS")
#write.csv( combined_df , "combined_df.csv")
saveRDS( combined_df , "combined_df.rds")
#' Save to csv for viewing in Excel and exporting to presentation (trim columns and rows to reduce size)
#' order by all waves with sufficient data and save top 1000 leading indicators
combined_df_order_perc_rank_w3_4_5_6_7 = combined_df[ with( combined_df , order( -combined_df$w3_4_5_6_7_perc_rank_mean ) ) , ]
plot( combined_df_order_perc_rank_w3_4_5_6_7$w3_4_5_6_7_perc_rank_mean )
#View( t( combined_df_order_perc_rank_w3_4_5_6_7[1:1000,] ) )
#View( combined_df_order_perc_rank_w3_4_5_6_7[ 1:1000 , c( 1:9 ,22,23,36,37,50,51,64,65,78,79,92,93,106,107, 150:167 ) ])
write.csv( combined_df_order_perc_rank_w3_4_5_6_7[ 1:5000 , c( 1:9 ,22,23,36,37,50,51,64,65,78,79,92,93,106,107, 150:167 ) ] , "combined_df_order_perc_rank_w3_4_5_6_7_top1000.csv")
#' order by genomic driven waves with sufficient data and save top 1000 leading indicators
combined_df_order_perc_rank_w3_4_7 = combined_df[ with( combined_df , order( -combined_df$w3_4_7_perc_rank_mean ) ) , ]
plot( combined_df_order_perc_rank_w3_4_7$w3_4_7_perc_rank_mean )
write.csv( combined_df_order_perc_rank_w3_4_7[ 1:5000 , c( 1:9 ,22,23,36,37,50,51,64,65,78,79,92,93,106,107, 150:167 ) ] , "combined_df_order_perc_rank_w3_4_7_top1000.csv")
#' order by non-genomic driven waves with sufficient data and save top 1000 leading indicators
combined_df_order_perc_rank_w5_6 = combined_df[ with( combined_df , order( -combined_df$w5_6_perc_rank_mean ) ) , ]
plot( combined_df_order_perc_rank_w5_6$w5_6_perc_rank_mean )
write.csv( combined_df_order_perc_rank_w5_6[ 1:5000 , c( 1:9 ,22,23,36,37,50,51,64,65,78,79,92,93,106,107, 150:167 ) ] , "combined_df_order_perc_rank_w5_6_top1000.csv")


#'** actually probably better not to filter as will skew the percentile rankings**
combined_df_filter_normMCC_07$w_all_ex_2_8_normMCC = round( combined_df_filter_normMCC_07$w_all_ex_2_8_normMCC , 3 )
combined_df_filter_normMCC_07$w_all_ex_2_8_ROC_AUC = round( combined_df_filter_normMCC_07$w_all_ex_2_8_ROC_AUC , 3 )
combined_df_filter_normMCC_07$min_normMCC_w3_to_w7 = round( combined_df_filter_normMCC_07$min_normMCC_w3_to_w7 , 3 )
combined_df_filter_normMCC_07$geom_mean_normMCC_w3_to_w7  = round( combined_df_filter_normMCC_07$geom_mean_normMCC_w3_to_w7  , 3 )

combined_df_filter_normMCC_07_order_multi = combined_df_filter_normMCC_07[ with( combined_df_filter_normMCC_07 
                                                               , order(   -w_all_ex_2_8_normMCC
                                                                        , -w_all_ex_2_8_ROC_AUC
                                                                        , -min_normMCC_w3_to_w7
                                                                        , -geom_mean_normMCC_w3_to_w7 
                                                                        ) 
                                                               ) 
                                                         , ]

View( t( combined_df_filter_normMCC_07_order_multi ))
plot( combined_df_filter_normMCC_07_order_multi$w_all_ex_2_8_normMCC , combined_df_filter_normMCC_07_order_multi$min_normMCC_w3_to_w7)
max( combined_df_filter_normMCC_07_order_multi$min_normMCC_w3_to_w7 )
View( t( subset( combined_df_filter_normMCC_07_order_multi 
             , combined_df_filter_normMCC_07_order_multi$min_normMCC_w3_to_w7 == max( combined_df_filter_normMCC_07_order_multi$min_normMCC_w3_to_w7 ) )))

#' Select based on highest minimum normalised MCC across waves 3 to 7
message( "Highest minimum normalised MCC = " , round( max( combined_df$min_normMCC_w3_to_w7 ) ,3 ) )
print( t( subset( combined_df 
                 , combined_df$min_normMCC_w3_to_w7 == max( combined_df$min_normMCC_w3_to_w7 ) ))[1:9,])
#' Select based on highest geo normalised MCC across waves 3 to 7
message( "Highest geometric mean normalised MCC = " , round( max( combined_df$geom_mean_normMCC_w3_to_w7 ),3) )
print( t( subset( combined_df 
                 , combined_df$geom_mean_normMCC_w3_to_w7 == max( combined_df$geom_mean_normMCC_w3_to_w7 ) ))[1:9,])
#' Select based on highest arithmetic normalised MCC across waves 3 to 7
message( "Highest arithmetic mean normalised MCC = ", round(max( combined_df$arith_mean_normMCC_w3_to_w7 ),3) )
print( t( subset( combined_df 
                 , combined_df$arith_mean_normMCC_w3_to_w7 == max( combined_df$arith_mean_normMCC_w3_to_w7 ) ))[1:9,])
#' Select based on highest arithmetic normalised MCC across waves 3 to 7
message( "Highest ROC AUC = ", round( max( combined_df$w_all_ex_2_8_ROC_AUC ),3) )
print( t( subset( combined_df 
                 , combined_df$w_all_ex_2_8_ROC_AUC == max( combined_df$w_all_ex_2_8_ROC_AUC ) ))[1:9,])
#' Select based on highest normalised MCC across waves 3 to 7
message( "Highest normalised MCC = " , round( max( combined_df$w_all_ex_2_8_normMCC ), 3) )
print( t( subset( combined_df 
                 , combined_df$w_all_ex_2_8_normMCC == max( combined_df$w_all_ex_2_8_normMCC ) ))[1:9,]
       )


#' Show top x leading indicators when ordered by different parameters
x = 20
order_normMCC = combined_df[ with( combined_df
                         , order(   -w_all_ex_2_8_normMCC
                                  #, -w_all_ex_2_8_ROC_AUC
                                  #, -min_normMCC_w3_to_w7
                                  #, -geom_mean_normMCC_w3_to_w7 
                                    ) ), ][ 1:x , c( 1:9 , 148:152 ) ]
order_normMCC$code = paste( order_normMCC$time_window_anchor , order_normMCC$type , order_normMCC$leading_indicator , order_normMCC$min_age, order_normMCC$max_age, order_normMCC$min_desc, order_normMCC$parent_sub_lgr_threshold, order_normMCC$pval_threshold, order_normMCC$ews_threshold)
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/TFPS vs non-TFPS")
write.csv( order_normMCC , "top20_by_normMCC")
order_min = combined_df[ with( combined_df
                         , order(   -min_normMCC_w3_to_w7
                                    #, -w_all_ex_2_8_normMCC
                                    #, -w_all_ex_2_8_ROC_AUC
                                    #, -geom_mean_normMCC_w3_to_w7 
                         ) ), ][1:x,c(1:9,148:152)]
order_min$code = paste( order_min$time_window_anchor , order_min$type , order_min$leading_indicator , order_min$min_age, order_min$max_age, order_min$min_desc, order_min$parent_sub_lgr_threshold, order_min$pval_threshold, order_min$ews_threshold)

order_geo_mean = combined_df[ with( combined_df
                               , order(   -geom_mean_normMCC_w3_to_w7
                                          #, -w_all_ex_2_8_normMCC
                                          #, -w_all_ex_2_8_ROC_AUC
                                          #,  ,-min_normMCC_w3_to_w7
                               ) ), ][1:x,c(1:9,148:152)]
order_geo_mean$code = paste( order_geo_mean$time_window_anchor , order_geo_mean$type , order_geo_mean$leading_indicator , order_geo_mean$min_age, order_geo_mean$max_age, order_geo_mean$min_desc, order_geo_mean$parent_sub_lgr_threshold, order_geo_mean$pval_threshold, order_geo_mean$ews_threshold )

order_ROC_AUC = combined_df[ with( combined_df
                                    , order(   -w_all_ex_2_8_ROC_AUC
                                            ) ), ][1:x,c(1:9,148:152)]
order_ROC_AUC$code = paste( order_ROC_AUC$time_window_anchor , order_ROC_AUC$type , order_ROC_AUC$leading_indicator , order_ROC_AUC$min_age, order_ROC_AUC$max_age, order_ROC_AUC$min_desc, order_ROC_AUC$parent_sub_lgr_threshold, order_ROC_AUC$pval_threshold, order_ROC_AUC$ews_threshold)

order_search = rbind( order_normMCC$code , order_min$code , order_geo_mean$code , order_ROC_AUC$code )
View(table(order_search))

View( t(subset( combined_df , combined_df$time_window_anchor == "Rt+5" &
                  combined_df$type == "Genomic" &
                  combined_df$leading_indicator == "tfps_lgr_mean" &
                  combined_df$min_age == 28 &
                  combined_df$max_age == 84 &
                  combined_df$min_desc == "perc" &
                  combined_df$parent_sub_lgr_threshold == "075" &
                  combined_df$pval_threshold == "005" &
                  combined_df$ews_threshold == 0.5
))
)


plot( combined_df_filter_normMCC_07_order_multi$w_all_ex_2_8_normMCC , combined_df_filter_normMCC_07_order_multi$min_normMCC_w3_to_w7)
max( combined_df_filter_normMCC_07_order_multi$min_normMCC_w3_to_w7 )


#' Look at stats for best genomic leading indicators using TFPS method
View( t(subset( combined_df , combined_df$leading_indicator == "tfps_pango_lineage_max_lgr" &
                combined_df$min_age == 7 &
                combined_df$max_age == 56 &
                combined_df$min_desc == 20 &
                combined_df$parent_sub_lgr_threshold == "085" &
                combined_df$pval_threshold == "001" &
                combined_df$ews_threshold == 0
              )
      ))
View( t(subset( combined_df , combined_df$leading_indicator == "tfps_lgr_simple_mean" &
                combined_df$min_age == 14 &
                combined_df$max_age == 56 &
                combined_df$min_desc == 20 &
                combined_df$parent_sub_lgr_threshold == "085" &
                combined_df$pval_threshold == "no" &
                combined_df$ews_threshold == 0
))
)
rank1 = subset( combined_df_order_perc_rank_w3_4_5_6_7 , combined_df_order_perc_rank_w3_4_5_6_7$leading_indicator == "tfps_pango_lineage_max_lgr" &
                  combined_df_order_perc_rank_w3_4_5_6_7$min_age == 7 &
                  combined_df_order_perc_rank_w3_4_5_6_7$max_age == 56 &
                  combined_df_order_perc_rank_w3_4_5_6_7$min_desc == 20 &
                  combined_df_order_perc_rank_w3_4_5_6_7$parent_sub_lgr_threshold == "085" &
                  combined_df_order_perc_rank_w3_4_5_6_7$pval_threshold == "001" &
                  combined_df_order_perc_rank_w3_4_5_6_7$ews_threshold == 0
              )
View(t(rank1))
rank2 = subset( combined_df_order_perc_rank_w3_4_5_6_7 , combined_df_order_perc_rank_w3_4_5_6_7$leading_indicator == "tfps_lgr_simple_mean" &
                  combined_df_order_perc_rank_w3_4_5_6_7$min_age == 14 &
                  combined_df_order_perc_rank_w3_4_5_6_7$max_age == 56 &
                  combined_df_order_perc_rank_w3_4_5_6_7$min_desc == 20 &
                  combined_df_order_perc_rank_w3_4_5_6_7$parent_sub_lgr_threshold == "085" &
                  combined_df_order_perc_rank_w3_4_5_6_7$pval_threshold == "no" &
                  combined_df_order_perc_rank_w3_4_5_6_7$ews_threshold == 0
              )
View(t(rank2))
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/TFPS vs non-TFPS")
write.csv( rank1[ , c( 1:9 ,22,23,36,37,50,51,64,65,78,79,92,93,106,107, 150:167 ) ] , "rank1_genomic_ranking.csv")
write.csv( rank2[ , c( 1:9 ,22,23,36,37,50,51,64,65,78,79,92,93,106,107, 150:167 ) ] , "rank2_genomic_ranking.csv")



#' Genomic TFPS. start date time anchor, -30<t<+10d
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/ROC_multi_stat_HPC_combined/filtered/" )
saveRDS( data_filtered_AUC_07_order , "ROC_stats_sd_tfps_10d_filtered_065.rds")
write.csv( data_filtered_AUC_07_order , "ROC_stats_sd_tfps_10d_filtered_065.csv")
#' Genomic TFPS. start date time anchor, -30<t<+5d
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/analysis/ROC_multi_stat/ROC_multi_stat_HPC_combined/filtered/" )
saveRDS( data_filtered_AUC_07_order , "ROC_stats_sd_tfps_5d_filtered_065.rds")
write.csv( data_filtered_AUC_07_order , "ROC_stats_sd_tfps_5d_filtered_065.csv")
#' Non-genomic. Rt critical transition time anchor, -30<t<+5d
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/filtered/" )
saveRDS( data_filtered_AUC_07_order , "ROC_stats_rt_non_tfps_5d_filtered_065.rds")
write.csv( data_filtered_AUC_07_order , "ROC_stats_rt_non_tfps_5d_filtered_065.csv")
#' Non-genomic. Rt critical transition time anchor, -30<t<+10d
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/filtered/" )
saveRDS( data_filtered_AUC_07_order , "ROC_stats_rt_non_tfps_10d_filtered_065.rds")
write.csv( data_filtered_AUC_07_order , "ROC_stats_rt_non_tfps_10d_filtered_065.csv")
#' Non-genomic. wave start date time anchor, -30<t<+5d
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/filtered/" )
saveRDS( data_filtered_AUC_07_order , "ROC_stats_sd_non_tfps_5d_filtered_065.rds")
write.csv( data_filtered_AUC_07_order , "ROC_stats_sd_non_tfps_5d_filtered_065.csv")
#' Non-genomic. wave start date time anchor, -30<t<+10d
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/filtered/" )
saveRDS( data_filtered_AUC_07_order , "ROC_stats_sd_non_tfps_10d_filtered_065.rds")
write.csv( data_filtered_AUC_07_order , "ROC_stats_sd_non_tfps_10d_filtered_065.csv")

#' Look at all leading indicator parameter sets for Google mobility workplaces and 7d rolling average
combined_df_order_perc_rank_w3_4_5_6_7_workplaces = subset( combined_df_order_perc_rank_w3_4_5_6_7 
              ,  leading_indicator == "workplaces_percent_change_from_baseline" |
                 leading_indicator == "workplaces_percent_change_from_baseline_7d_mean" 
              )
View( combined_df_order_perc_rank_w3_4_5_6_7_workplaces )
temp = subset( combined_df_order_perc_rank_w3_4_5_6_7_workplaces, 
               abs( ews_threshold - 0.30) < 0.001 ) 
View( temp )

#' Look at all leading indicator parameter sets where normalised MCC is not defined
#' wave 3
temp = subset( non_tfps_rt_5d 
              , w3_normMCC == "NaN" )
test = paste0( temp$w3_TP,"_", temp$w3_FP,"_", temp$w3_TN,"_", temp$w3_FN )
unique( test )  
temp = subset( non_tfps_rt_5d 
               , w4_normMCC == "NaN" )
test = paste0( temp$w4_TP,"_", temp$w4_FP,"_", temp$w4_TN,"_", temp$w4_FN )
unique( test )  

#**OLD - USED FOR AUC ONLY RESULTS**
#' Plotting
AUC_df = AUC_rt_tfps_combined_5d #AUC_rt_tfps_combined_10d , AUC_sd_tfps_combined_5d  AUC_sd_tfps_combined_10d

plot( subset( AUC_df , AUC_df$new_col_li_rt == unique( AUC_df$new_col_li_rt )[ 1 ] )[ , 16 ]
      , ylim=c(-1,1)
      , ylab = "AUC (x=FPR,y=TPR)"
      , xlab = "Different TFP scanner parameters and cluster filters"
      , col=rainbow(19)[1]
      ,cex.lab = 1.5
      ,cex.axis = 1.5)
abline( h = 0 ) ; abline( h = c(-0.7,0.7) , lty = 2 )
#' plot remaining leading indicators
for(i in 2:19){
  points( subset( AUC_df , AUC_df$new_col_li_rt == unique( AUC_df$new_col_li_rt )[ i ] )[ , 16 ]
          , col=rainbow(19)[i]
  )
}


#' Mean plot
AUC_df_mean = mean( as.numeric( subset( AUC_df , AUC_df$new_col_li_rt == unique( AUC_df$new_col_li_rt )[ 1 ] )[ , 16 ] ) , na.rm=TRUE) 
plot( seq(1,nrow(AUC_df),1)
      , replicate( nrow(AUC_df) , AUC_df_mean )
      , ylim=c(0,1)
      , ylab = "AUC (x=FPR,y=TPR)"
      , xlab = "Different TFP scanner parameters and cluster filters"
      , col=rainbow(19)[1]
      ,cex.lab = 1.5
      ,cex.axis = 1.5
      ,typ="l")
abline( h = 0 ) ; abline( h = c(-0.7,0.7) , lty = 2 )
#' plot remaining leading indicators
for(i in 2:19){
  AUC_df_mean = mean( as.numeric( subset( AUC_df , AUC_df$new_col_li_rt == unique( AUC_df$new_col_li_rt )[ i ] )[ , 16 ] ) , na.rm=TRUE) 
  message(i," ",round(AUC_df_mean,2))
  points( x = seq( 1 , nrow( AUC_df ) , 1 )
          , y = replicate( nrow( AUC_df ) , AUC_df_mean )
          , col=rainbow(19)[i]
          , typ="l"
  )
}


legend("topright" , legend=c( unique(AUC_df$new_col_li_rt ) ) , col = c(rainbow(19)) ,pch=replicate(19,16) )

max_auc = max(AUC_df[,16]) ; max_index = which( AUC_df[,16] == max_auc) ; AUC_df[max_index,]
min(AUC_df[,16])
rank_level = 0.6
auc_sub = subset( AUC_df , AUC_df$wave_3 >= rank_level &
                    AUC_df$wave_4 >= rank_level &
                    AUC_df$wave_5 >= rank_level &
                    AUC_df$wave_6 >= rank_level &
                    AUC_df$wave_7 >= rank_level &
                    AUC_df$all_exc_B.1.177_BA.1 >= rank_level
)
View(auc_sub)
