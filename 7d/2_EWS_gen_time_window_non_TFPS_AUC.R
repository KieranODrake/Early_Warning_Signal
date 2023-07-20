#' Generate EWS from non-TFPS data using different method from that applied to main analysis of 
#' TFPS derived leading indicators
#' Purpose: Generate early warning signals (EWS) using leading indicator data (non-TFPS, for TFPS EWS see 'EWS_calc_threshold_lgr_th_XXX_pval_filter_XXX_HPC_array.R')
#' A range of leading indicator thresholds are used for generating EWS. 
#' Requires data to be loaded by running this script first 
#' '1_EWS_gen_time_window_non_TFPS_processing.R'
#' Calculates the area under the curve (AUC) of the ROC for each leading indicator type.
#' See '2_EWS_gen_time_window_non_TFPS_multi_stat.R' for calculation of other ROC stats which may be more appropriate
#' See Chicco and Jurman BioData Mining (2023) 16:4  https://doi.org/10.1186/s13040-023-00322-4
#' Author: Kieran Drake
#' Date: March 2023



#' 6 - Can then calculate AUC and various other ROC stats (https://en.wikipedia.org/wiki/Receiver_operating_characteristic)
#'     for each wave and for whole period (average or total TP, FP, TN, FN across all waves).

library( magrittr )
library( lubridate )
#library( ggplot2 )
library( stringr )
library( zoo )
library( data.table )
library( sqldf )
library( zoo )




#' Define area under curve calculation (AUC) function using trapezoidal rule of integration
#' https://en.wikipedia.org/wiki/Trapezoidal_rule and https://stackoverflow.com/questions/7358738/how-to-calculate-integral-with-r
AUC <- function(x, y){
  sum(diff(x)*rollmean(y,2))
}

#ews_class_total = table( li_z_ews_rt_classification_array[ c( date_index_total ) , 1 , 201 ] )

#' Create dataframes of AUC values
AUC_rt_df = data.frame( "wave" = c("wave_2","wave_3","wave_4","wave_5","wave_6","wave_7","wave_8","all","all_exc_B.1.177","all_exc_B.1.177_BA.1") )
AUC_sd_df = data.frame( "wave" = c("wave_2","wave_3","wave_4","wave_5","wave_6","wave_7","wave_8","all","all_exc_B.1.177","all_exc_B.1.177_BA.1") )
TPR_rt = data.frame(replicate(dim( li_z_ews_rt_classification_array )[ 3 ],NA)); FPR_rt = TPR_rt; TPR_sd = TPR_rt; FPR_sd = TPR_rt ; MCC_rt = TPR_rt ; MCC_sd = TPR_sd #' Initialise
date_index_list = list(date_index_w2,date_index_w3,date_index_w4,date_index_w5,date_index_w6,date_index_w7,date_index_w8
                      ,date_index_total,date_index_total_ex_B.1.177,date_index_total_ex_B.1.177_BA.1 )
names(date_index_list) <- c("date_index_w2","date_index_w3","date_index_w4","date_index_w5","date_index_w6","date_index_w7"
                            ,"date_index_w8","date_index_total","date_index_total_ex_B.1.177","date_index_total_ex_B.1.177_BA.1" )
par(mfrow=c(2,1)) #' set up plot area

###** For testing **
#li=20
#date_index = date_index_list[[7]]
#' Then jump to ews_th for loop
#####################

#' Calculation
#' loop through leading indicators
for( li in 1 : dim( li_z_ews_rt_classification_array )[2] ){
  row_counter = 0
  #' loop through wave date ranges inc all waves, all waves exc. B.1.177 and all waves exc B.1.177 & BA.2
  #' and calculate ROC stats and AUC for each wave etc.
  for( date_index in date_index_list ){
    row_counter = row_counter + 1
    wave_name = c("wave 2 B.1.177","wave 3 Alpha","wave 4 Delta(1)","wave 5 Delta(2)","wave 6 Delta(3)","wave 7 BA.1","wave 8 BA.2","all waves","all waves exc. B.1.177","all waves exc. B.1.177 & BA.1")
    wave_name = wave_name[ row_counter - 1 ]
    for( ews_th in 1 : dim( li_z_ews_rt_classification_array )[ 3 ] ){
      #' Using Rt critical transition date
      ews_class_total_rt = table( li_z_ews_rt_classification_array[ date_index , li , ews_th ] )
      TP_rt = ifelse( is.na( ews_class_total_rt["TP"] ) , 0 , ews_class_total_rt["TP"] )
      FP_rt = ifelse( is.na( ews_class_total_rt["FP"] ) , 0 , ews_class_total_rt["FP"] )
      TN_rt = ifelse( is.na( ews_class_total_rt["TN"] ) , 0 , ews_class_total_rt["TN"] )
      FN_rt = ifelse( is.na( ews_class_total_rt["FN"] ) , 0 , ews_class_total_rt["FN"] )
      TPR_rt[ ews_th , 1 ] = ifelse( is.na( TP_rt / (TP_rt + FN_rt) ) , 0 , TP_rt / (TP_rt + FN_rt) )
      FPR_rt[ ews_th , 1 ] = ifelse( is.na( FP_rt / (FP_rt + TN_rt) ) , 0 , FP_rt / (FP_rt + TN_rt) )
      MCC_rt[ ews_th , 1 ] = ((TP_rt*TN_rt) - (FP_rt*FN_rt)) / sqrt( (TP_rt+FP_rt)*(TP_rt+FN_rt)*(TN_rt+FP_rt)*(TN_rt+FN_rt) )
      #' Using wave start(inflection) date based on hospitalisations
      ews_class_total_sd = table( li_z_ews_sd_classification_array[ date_index , li , ews_th ] )
      TP_sd = ifelse( is.na( ews_class_total_sd["TP"] ) , 0 , ews_class_total_sd["TP"] )
      FP_sd = ifelse( is.na( ews_class_total_sd["FP"] ) , 0 , ews_class_total_sd["FP"] )
      TN_sd = ifelse( is.na( ews_class_total_sd["TN"] ) , 0 , ews_class_total_sd["TN"] )
      FN_sd = ifelse( is.na( ews_class_total_sd["FN"] ) , 0 , ews_class_total_sd["FN"] )
      TPR_sd[ ews_th , 1 ] = ifelse( is.na( TP_sd / (TP_sd + FN_sd) ) , 0 , TP_sd / (TP_sd + FN_sd) )
      FPR_sd[ ews_th , 1 ] = ifelse( is.na( FP_sd / (FP_sd + TN_sd) ) , 0 , FP_sd / (FP_sd + TN_sd) ) 
      MCC_sd[ ews_th , 1 ] = ((TP_sd*TN_sd) - (FP_sd*FN_sd)) / sqrt( (TP_sd+FP_sd)*(TP_sd+FN_sd)*(TN_sd+FP_sd)*(TN_sd+FN_sd) )
      #sensitivity = TPR
      #specificity = TN / (TN + FP)
      #youdens_j = sensitivity + specificity -1
    }
    TPR_FPR_rt = data.frame("FPR"=FPR_rt , "TPR"=TPR_rt) ; names(TPR_FPR_rt)[c(1,2)] <- c("FPR","TPR") ; TPR_FPR_rt = TPR_FPR_rt[ with( TPR_FPR_rt , order( FPR , TPR ) ), ]
    TPR_FPR_sd = data.frame("FPR"=FPR_sd , "TPR"=TPR_sd) ; names(TPR_FPR_sd)[c(1,2)] <- c("FPR","TPR") ; TPR_FPR_sd = TPR_FPR_sd[ with( TPR_FPR_sd , order( FPR , TPR ) ), ]
    #names( table(li_z_ews_rt_classification_array[,1,201]  ))#[1]
    ROC_AUC_rt = AUC(x=TPR_FPR_rt$FPR, y=TPR_FPR_rt$TPR)
    ROC_AUC_sd = AUC(x=TPR_FPR_sd$FPR, y=TPR_FPR_sd$TPR)
    #' Plot ROC for Rt critical transition dates
    plot( TPR_FPR_rt$FPR , TPR_FPR_rt$TPR , typ="o" , col="blue" , main = paste0( colnames( li_z_df )[ li+1 ] , " | ", wave_name ) ) 
    mtext( paste0("AUC = ",round( ROC_AUC_rt , 2 )) , at = c( 0.5 ) )
    lines(seq(0,1,0.01),seq(0,1,0.01),typ="l")
    #' Plot ROC for wave start (inflection) dates
    plot( TPR_FPR_sd$FPR , TPR_FPR_sd$TPR , typ="o" , col="red" , main = paste0( colnames( li_z_df )[ li+1 ] , " | ", wave_name ) ) 
    mtext( paste0("AUC = ",round( ROC_AUC_sd , 2 )) , at = c( 0.5 ) )
    lines(seq(0,1,0.01),seq(0,1,0.01),typ="l")
    #' Add AUC values to the respective dataframes
    AUC_rt_df[ row_counter , li+1 ] = ROC_AUC_rt
    AUC_sd_df[ row_counter , li+1 ] = ROC_AUC_sd
  }
  #' Add leading indicator names to columns 
  names(AUC_rt_df)[li+1] <- names(li_z_df)[li+1]
  names(AUC_sd_df)[li+1] <- names(li_z_df)[li+1]
}

#' Plot example of ROC curves
#' 
par(mfrow=c(1,1))
plot( 0 , 0 , xlim=c(0,1) , ylim=c(0,1)
      ,main = paste0( "Time window around Rt critical transition")# | all waves exc. B.1.177 & BA.2" )
      ,xlab="False positivity rate (FPR) (1 - specificity) ",ylab="True positivity rate (TPR) (Sensitivity)"
      ) 
date_index = date_index_list[[10]]
li=max_index
#' loop through leading indicators
for( li in 1 : dim( li_z_ews_rt_classification_array )[2] ){
  row_counter = 0
  #' loop through wave date ranges inc all waves, all waves exc. B.1.177 and all waves exc B.1.177 & BA.2
  #' and calculate ROC stats and AUC for each wave etc.
  for( date_index in date_index_list ){
    row_counter = row_counter + 1
    wave_name = c("wave 2 B.1.177","wave 3 Alpha","wave 4 Delta(1)","wave 5 Delta(2)","wave 6 Delta(3)","wave 7 BA.1","wave 8 BA.2","all waves","all waves exc. B.1.177","all waves exc. B.1.177 & BA.1")
    wave_name = wave_name[ row_counter - 1 ]
    for( ews_th in 1 : dim( li_z_ews_rt_classification_array )[ 3 ] ){
      #' Using Rt critical transition date
      ews_class_total_rt = table( li_z_ews_rt_classification_array[ date_index , li , ews_th ] )
      TP_rt = ifelse( is.na( ews_class_total_rt["TP"] ) , 0 , ews_class_total_rt["TP"] )
      FP_rt = ifelse( is.na( ews_class_total_rt["FP"] ) , 0 , ews_class_total_rt["FP"] )
      TN_rt = ifelse( is.na( ews_class_total_rt["TN"] ) , 0 , ews_class_total_rt["TN"] )
      FN_rt = ifelse( is.na( ews_class_total_rt["FN"] ) , 0 , ews_class_total_rt["FN"] )
      TPR_rt[ ews_th , 1 ] = ifelse( is.na( TP_rt / (TP_rt + FN_rt) ) , 0 , TP_rt / (TP_rt + FN_rt) )
      FPR_rt[ ews_th , 1 ] = ifelse( is.na( FP_rt / (FP_rt + TN_rt) ) , 0 , FP_rt / (FP_rt + TN_rt) )
      #' Using wave start(inflection) date based on hospitalisations
      ews_class_total_sd = table( li_z_ews_sd_classification_array[ date_index , li , ews_th ] )
      TP_sd = ifelse( is.na( ews_class_total_sd["TP"] ) , 0 , ews_class_total_sd["TP"] )
      FP_sd = ifelse( is.na( ews_class_total_sd["FP"] ) , 0 , ews_class_total_sd["FP"] )
      TN_sd = ifelse( is.na( ews_class_total_sd["TN"] ) , 0 , ews_class_total_sd["TN"] )
      FN_sd = ifelse( is.na( ews_class_total_sd["FN"] ) , 0 , ews_class_total_sd["FN"] )
      TPR_sd[ ews_th , 1 ] = ifelse( is.na( TP_sd / (TP_sd + FN_sd) ) , 0 , TP_sd / (TP_sd + FN_sd) )
      FPR_sd[ ews_th , 1 ] = ifelse( is.na( FP_sd / (FP_sd + TN_sd) ) , 0 , FP_sd / (FP_sd + TN_sd) ) 
      #sensitivity = TPR
      #specificity = TN / (TN + FP)
      #youdens_j = sensitivity + specificity -1
    }
    TPR_FPR_rt = data.frame("FPR"=FPR_rt , "TPR"=TPR_rt) ; names(TPR_FPR_rt)[c(1,2)] <- c("FPR","TPR") ; TPR_FPR_rt = TPR_FPR_rt[ with( TPR_FPR_rt , order( FPR , TPR ) ), ]
    TPR_FPR_sd = data.frame("FPR"=FPR_sd , "TPR"=TPR_sd) ; names(TPR_FPR_sd)[c(1,2)] <- c("FPR","TPR") ; TPR_FPR_sd = TPR_FPR_sd[ with( TPR_FPR_sd , order( FPR , TPR ) ), ]
    #names( table(li_z_ews_rt_classification_array[,1,201]  ))#[1]
    ROC_AUC_rt = AUC(x=TPR_FPR_rt$FPR, y=TPR_FPR_rt$TPR)
    ROC_AUC_sd = AUC(x=TPR_FPR_sd$FPR, y=TPR_FPR_sd$TPR)
    #' Plot ROC for Rt critical transition dates
    #plot( TPR_FPR_rt$FPR , TPR_FPR_rt$TPR , typ="o" , col="blue" , main = paste0( colnames( li_z_df )[ li+1 ] , " | ", wave_name ) ) 
    #lines( TPR_FPR_rt$FPR , TPR_FPR_rt$TPR , typ="l" , col=rainbow(71)[li] ) 
    lines( TPR_FPR_rt$FPR , TPR_FPR_rt$TPR , typ="o" , col=rainbow(10)[row_counter] ) 
    #mtext( paste0("AUC = ",round( ROC_AUC_rt , 2 )) , at = c( 0.5 ) )
    lines(seq(0,1,0.01),seq(0,1,0.01),typ="l")
    #' Plot ROC for wave start (inflection) dates
    #plot( TPR_FPR_sd$FPR , TPR_FPR_sd$TPR , typ="o" , col="red" , main = paste0( colnames( li_z_df )[ li+1 ] , " | ", wave_name ) ) 
    #mtext( paste0("AUC = ",round( ROC_AUC_sd , 2 )) , at = c( 0.5 ) )
    #lines(seq(0,1,0.01),seq(0,1,0.01),typ="l")
  }
}
legend()


#' Further look at best leading indicator from AUC point of view
plot(trim_comix_na_approx[,1],trim_comix_na_approx[,6],typ="l",ylim=c(0,5),ylab="mean daily contacts and new hospital admissions/200",xlab="")
lines(dat_df$date,dat_df$cases/1000,col="red")
lines(trim_comix_na_approx[,1],trim_comix_na_approx[,4],typ="l",col="blue")
lines(trim_comix_na_approx[,1]+40,trim_comix_na_approx[,4],typ="l",col="blue")

#' plot AUCs
par(mfrow=c(2,1)) #' set up plot area
wave_to_plot = 10 # individual waves = 1 2 3 4 5 6 7, total = 8, total ex B.1.177 = 9
hosp_ix = c(seq(2,5,1)) ; pos_rate_ix = c(seq(6,7,1)) ; Ct_mean_ix = c(seq(8,26,1)) ; Ct_median_ix = c(seq(27,45,1)) ; CoMix_ix = c(seq(46,58,1)) ; Google_mob_ix = c(seq(59,65,1)) ; Google_mob_7d_avg_ix = c(seq(66,72,1))
#' Using Rt critical transition as anchor for time window
plot(seq(2,72,1),AUC_rt_df[wave_to_plot,2:72],ylim=c(0,1),xlab="leading indicator column number",ylab="AUC (x=FPR,y=TPR)",main="Time window around Rt critical transition")
points( hosp_ix , AUC_rt_df[ wave_to_plot , hosp_ix ] , col=rainbow(7)[1] , pch = 16) ; points( pos_rate_ix , AUC_rt_df[ wave_to_plot , pos_rate_ix ] , col=rainbow(7)[2] , pch = 16 ) ; points( Ct_mean_ix , AUC_rt_df[ wave_to_plot , Ct_mean_ix ] , col=rainbow(7)[3] , pch = 16 ); points( Ct_median_ix , AUC_rt_df[ wave_to_plot , Ct_median_ix ] , col=rainbow(7)[4] , pch = 16 ); points( CoMix_ix , AUC_rt_df[ wave_to_plot , CoMix_ix ] , col=rainbow(7)[5] , pch = 16 ); points( Google_mob_ix , AUC_rt_df[ wave_to_plot , Google_mob_ix ] , col=rainbow(7)[6] , pch = 16 ); points( Google_mob_7d_avg_ix , AUC_rt_df[ wave_to_plot , Google_mob_7d_avg_ix ] , col=rainbow(7)[7] , pch = 16 )
abline(h=0.5) ; abline(h=0.7,lty=2) ; abline(h=0.3,lty=2)
max_auc = max(AUC_rt_df[wave_to_plot,2:72]) ; max_index = which( AUC_rt_df[wave_to_plot,] == max_auc) ; colnames(AUC_rt_df)[max_index]
min(AUC_rt_df[wave_to_plot,2:72])
#' Using wave start/inflection date as anchor for time window
plot(seq(2,72,1),AUC_sd_df[wave_to_plot,2:72],ylim=c(0,1),xlab="leading indicator column number",ylab="AUC (x=FPR,y=TPR)",main="Time window around wave start/inflection date")
points( hosp_ix , AUC_sd_df[ wave_to_plot , hosp_ix ] , col=rainbow(7)[1] , pch = 16) ; points( pos_rate_ix , AUC_sd_df[ wave_to_plot , pos_rate_ix ] , col=rainbow(7)[2] , pch = 16 ) ; points( Ct_mean_ix , AUC_sd_df[ wave_to_plot , Ct_mean_ix ] , col=rainbow(7)[3] , pch = 16 ); points( Ct_median_ix , AUC_sd_df[ wave_to_plot , Ct_median_ix ] , col=rainbow(7)[4] , pch = 16 ); points( CoMix_ix , AUC_sd_df[ wave_to_plot , CoMix_ix ] , col=rainbow(7)[5] , pch = 16 ); points( Google_mob_ix , AUC_sd_df[ wave_to_plot , Google_mob_ix ] , col=rainbow(7)[6] , pch = 16 ); points( Google_mob_7d_avg_ix , AUC_sd_df[ wave_to_plot , Google_mob_7d_avg_ix ] , col=rainbow(7)[7] , pch = 16 )
abline(h=0.5) ; abline(h=0.7,lty=2) ; abline(h=0.3,lty=2)

legend("topright"
       ,legend=c("Hospitalisations","Positivity rate","PCR Ct - mean","PCR Ct - median","CoMix","Google mobility","Google mobility: 7d mean")
       ,col=c(rainbow(7)[1],rainbow(7)[2],rainbow(7)[3],rainbow(7)[4],rainbow(7)[5],rainbow(7)[6],rainbow(7)[7])
       ,pch=c(16,16,16,16,16,16,16)
       )

#' Change format of dataframes
AUC_sd_df_t = t(AUC_sd_df) ; colnames(AUC_sd_df_t) = AUC_sd_df_t[1,] ; AUC_sd_df_t = AUC_sd_df_t[-1,] ; AUC_sd_df_t = as.data.frame(AUC_sd_df_t)
AUC_rt_df_t = t(AUC_rt_df) ; colnames(AUC_rt_df_t) = AUC_rt_df_t[1,] ; AUC_rt_df_t = AUC_rt_df_t[-1,] ; AUC_rt_df_t = as.data.frame(AUC_rt_df_t)

#' Save dataframes
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/AUC/" )
saveRDS( AUC_sd_df_t , "AUC_sd_non_tfps_5d.rds")
saveRDS( AUC_rt_df_t , "AUC_rt_non_tfps_5d.rds")

AUC_sd_non_tfps_10d

#' 6 - Can then calculate AUC and various other ROC stats (https://en.wikipedia.org/wiki/Receiver_operating_characteristic)
#'     for each wave and for whole period (average or total TP, FP, TN, FN across all waves).
