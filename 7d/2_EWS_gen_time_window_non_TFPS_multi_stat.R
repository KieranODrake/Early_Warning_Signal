#' Purpose: Calculate receiver operating characteristic (ROC) stats 
#' (https://en.wikipedia.org/wiki/Receiver_operating_characteristic)
#' for early warning signals (EWS) generated from non-TFPS data 
#' (non-TFPS, for TFPS EWS see 'EWS_calc_threshold_lgr_th_XXX_pval_filter_XXX_HPC_array.R')
#' using different method from that applied to main analysis of TFPS derived leading indicators.
#' A range of leading indicator thresholds are used for generating EWS. 
#' Requires data to be loaded by running this script first 
#' '1_EWS_gen_time_window_non_TFPS_processing.R'
#' Calculates various binary classifier stats for each leading indicator type.
#' See '2_EWS_gen_time_window_non_TFPS_AUC.R' for calculation of ROC AUC only
#' But note Chicco & Jurman BioData Mining (2023) 16:4  https://doi.org/10.1186/s13040-023-00322-4
#' Author: Kieran Drake
#' Date: April 2023

library( magrittr )
library( lubridate )
library( stringr )
library( zoo )
library( data.table )
library( sqldf )

#' Define area under curve calculation (AUC) function using trapezoidal rule of integration
#' https://en.wikipedia.org/wiki/Trapezoidal_rule and https://stackoverflow.com/questions/7358738/how-to-calculate-integral-with-r
AUC <- function(x, y){
  sum(diff(x)*rollmean(y,2))
}

#' Calculation
#' 1 - Loop through the two different ways of calculating the time window (anchoring to Rt critical transition date or hospitalisation wave inflection date)
#' 2 - Loop through date ranges (individual waves, inc all waves, all waves exc. B.1.177 and all waves exc B.1.177 & BA.2)
#' 3 - Loop through leading indicators
#' 4 - Calculate ROC stats and AUC for each wave etc.
 
for( arr in 1:2 ){
  message(arr)
  #' initialise output data.frame template
  ROC_stat_types = c("TP","FP","TN","FN","TPR","FPR","TNR","PPV","NPV","F1","FMI","MCC","normMCC","ROC_AUC")
  output = data.frame( matrix( NA , nrow = ( ( ncol( li_z_df ) - 1 ) * dim( li_z_ews_rt_classification_array )[ 3 ]) #' <number of leading indicators> * <number of ews threshold values>
                               , ncol = 2 + ( length( ROC_stat_types ) * length( date_index_list ) ) #' number of date ranges plus 2 columns for the leading indicator type and the ews threshold value
                              ) 
                      ) 
  #' Add column names to output template
  colnames(output)[1:2] = c("leading_indicator_type","ews_threshold")
  counter=0
  for(wave_no in 2:(length(date_index_list)+1)){
    print(wave_no)
    for(ROC_stat in ROC_stat_types){
      counter = counter + 1
      if       (wave_no<= 8){ prefix <- c(paste0("w",wave_no,"_")) 
      } else if(wave_no== 9){ prefix <- c(paste0("w_all_"       ))
      } else if(wave_no==10){ prefix <- c(paste0("w_all_ex_2_"  )) 
      } else if(wave_no==11){ prefix <- c(paste0("w_all_ex_2_8_")) }
      new_col_name = paste0( prefix , ROC_stat )
      message(new_col_name)
      colnames(output)[counter+2] = new_col_name
    }
  }
  
  date_range_counter = 0
  for( date_index in date_index_list ){
    date_range_counter = date_range_counter + 1
    
    for( li in 1 : dim( li_z_ews_rt_classification_array )[ 2 ] ){
      #' Select between Rt critical transition date and wave start date
      if       (arr == 1){ temp_array = li_z_ews_rt_classification_array[ date_index , li ,  ] #' Using Rt critical transition as anchor point for time window
      } else if(arr == 2){ temp_array = li_z_ews_sd_classification_array[ date_index , li ,  ] }#' Using hospitalisation wave start date as anchor point for time window
      message( "Date range: ", date_range_counter , ". Leading indicator: ", li )
      
      #' Analyse EWS across all EWS thresholds for the particular date range (date_index) and leading indicator (li)
      temp_val_list = list()
      for(i in 1:dim( temp_array )[ 2 ]){
        temp_val = list( temp_array[  , i ] )
        temp_val_list = append( temp_val_list , temp_val )
      }
      ews_class_total = lapply( temp_val_list , table ) #( li_z_ews_rt_classification_array[ date_index , li , ews_th ] )
      out = data.frame()
      temp_output <- lapply( seq_along( ews_class_total ) 
                           , function( x , i ){   TP = as.numeric( ifelse( is.na( x[[i]]["TP"] ) , 0 , x[[i]]["TP"] ) ) #' TP
                                                  FP = as.numeric( ifelse( is.na( x[[i]]["FP"] ) , 0 , x[[i]]["FP"] ) ) #' FP
                                                  TN = as.numeric( ifelse( is.na( x[[i]]["TN"] ) , 0 , x[[i]]["TN"] ) ) #' TN
                                                  FN = as.numeric( ifelse( is.na( x[[i]]["FN"] ) , 0 , x[[i]]["FN"] ) ) #' FN
                                                  TPR = TP / (TP + FN) #' TPR True Positive Rate aka Sensitivity, true positive rate, recall, hit rate
                                                  FPR = FP / (FP + TN) #' FPR False Positive Rate = 1 - Specificity
                                                  TNR = TN / (TN + FP) #' TNR = True Negative Rate aka Specificity, selectivity, true negative rate
                                                  PPV = TP / (TP + FP)  #' PPV aka Precision, positive predictive value
                                                  NPV = TN / (TN + FN) #' NPV = Negative predictive value
                                                  F1 = (2*TP) / (2*TP + FP + FN) #' F1
                                                  FMI = sqrt( PPV * TPR )  #' Fowlkes-Mallows index
                                                  #' Need to break down MCC calculation as otherwise results in NA with message 'NAs produced by integer overflow'
                                                  #a = as.numeric(TP+FP) ; b = as.numeric(TP+FN) ; c = as.numeric(TN+FP) ; d = as.numeric(TN+FN)
                                                  #a = as.integer(TP+FP) ; b = as.integer(TP+FN) ; c = as.integer(TN+FP) ; d = as.integer(TN+FN)
                                                  #e = a*b ; f = e*c ; g = f*d
                                                  MCC = ((TP*TN) - (FP*FN)) / sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) #' cpmbines the 4 basic rates of confusion matrix
                                                  #MCC = ((TP*TN) - (FP*FN)) / sqrt( a*b ) #' cpmbines the 4 basic rates of confusion matrix
                                                  normMCC = (MCC + 1) /2
                                                  #' Add ROC stats to dataframe
                                                  out[ 1, c(seq(1,13,1))] = c(TP,FP,TN,FN,TPR,FPR,TNR,PPV,NPV,F1,FMI,MCC,normMCC)
                                                  return(out)
                                                  } 
                           , x = ews_class_total )
      temp_output_df = data.frame()
      for( j in 1: length( temp_output ) ){
        temp_output_df = rbind( temp_output_df , as.data.frame( temp_output[[j]] ) )
      }
      
      #' Calculate area under curve (AUC) for TPR vs FPR
      #' Rows to include must be for a single leading indicator but across all EWS thresholds
      TPR_FPR = data.frame( "FPR" = temp_output_df[ ,  6 ] , "TPR" = temp_output_df[ ,  5 ] )
      names( TPR_FPR )[c(1,2)] <- c("FPR","TPR") ; TPR_FPR = TPR_FPR[ with( TPR_FPR , order( FPR , TPR ) ) , ]
      ROC_AUC = AUC( x = TPR_FPR$FPR, y = TPR_FPR$TPR )
      temp_output_df[ , length( ROC_stat_types ) ] = ROC_AUC
      colnames( temp_output_df ) = ROC_stat_types
      #' Add data for this date range and leading indicator to the main dataframe
      #' First need to determine location in the dataframe
      row_start = 1 + ( (li-1) * dim( li_z_ews_rt_classification_array )[ 3 ] )
      row_end   = dim( li_z_ews_rt_classification_array )[ 3 ] + ( (li-1) * dim( li_z_ews_rt_classification_array )[ 3 ] )
      col_start = 3 + ( ( date_range_counter -1 ) * length(ROC_stat_types) )
      col_end   = 16 + ( ( date_range_counter -1 ) * length(ROC_stat_types) )
      #' Check source and destination dimensions match
      dim( output[ row_start:row_end , col_start:col_end ] ) == dim( temp_output_df )
      #' Then add data to final dataframe
      output[ row_start:row_end , col_start:col_end ] = temp_output_df
      
      if(date_range_counter == 2){
        output[ row_start:row_end ,  1 ] = colnames(li_z_df)[li+1]
        output[ row_start:row_end ,  2 ] = seq(-5,5,0.05)#[ews_th]
      }
    } #' li for loop end
  } #' date_index for loop end
  
  #' Change name of output dataframe depending on which dataset was used
  if       (arr == 1){ new_name = "output_rt" #; assign( output , output_rt ) } #' Using Rt critical transition as anchor point for time window
  } else if(arr == 2){ new_name = "output_sd" }#; assign( output , output_sd ) }#' Using hospitalisation wave start date as anchor point for time window
  assign( new_name , output )

}

#' Save dataframes
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/" )
saveRDS( output_sd , "ROC_stats_sd_non_tfps_5d.rds")
saveRDS( output_rt , "ROC_stats_rt_non_tfps_5d.rds")
write.csv( output_sd , "ROC_stats_sd_non_tfps_5d.csv")
write.csv( output_rt , "ROC_stats_rt_non_tfps_5d.csv")
#saveRDS( output_sd , "ROC_stats_sd_non_tfps_10d.rds")
#saveRDS( output_rt , "ROC_stats_rt_non_tfps_10d.rds")
#write.csv( output_sd , "ROC_stats_sd_non_tfps_10d.csv")
#write.csv( output_rt , "ROC_stats_rt_non_tfps_10d.csv")

#' Sorting and filtering dataframes in search for 'best' non-TFPS leading indicator
#' 
output_sd = readRDS( "ROC_stats_sd_non_tfps_5d.rds")
output_rt = readRDS( "ROC_stats_rt_non_tfps_5d.rds")
output_sd = readRDS( "ROC_stats_sd_non_tfps_10d.rds")
output_rt = readRDS( "ROC_stats_rt_non_tfps_10d.rds")

ROC_stats_base = output_sd
ROC_stats_base = output_rt #' output_rt or output_sd

ROC_stats_filtered = subset( ROC_stats_base
                             , ROC_stats_base$w3_TP > 0 &
                               ROC_stats_base$w4_TP > 0 &
                               ROC_stats_base$w5_TP > 0 &
                               ROC_stats_base$w6_TP > 0 &
                               ROC_stats_base$w7_TP > 0
                            )
ROC_stats_filtered = ROC_stats_filtered[ with( ROC_stats_filtered , order( w_all_ex_2_8_normMCC , decreasing = TRUE) ) , ]
#' Look at what is changing when ROC AUC remains high but normalised MCC falls
filter_val = 0.7
ROC_stats_filtered = subset( ROC_stats_filtered 
                             , ROC_stats_filtered$w_all_ex_2_8_ROC_AUC > filter_val #&
                               #ROC_stats_filtered$w3_ROC_AUC > filter_val &
                               #ROC_stats_filtered$w4_ROC_AUC > filter_val &
                               #ROC_stats_filtered$w5_ROC_AUC > filter_val &
                               #ROC_stats_filtered$w6_ROC_AUC > filter_val &
                               #ROC_stats_filtered$w7_ROC_AUC > filter_val
                             )
plot(  x = ROC_stats_base$w_all_ex_2_8_normMCC 
       , y = ROC_stats_base$w_all_ex_2_8_ROC_AUC
       , xlab = "Normalised MCC (Alpha, Delta (1,2,3), BA.1)"
       , ylab = "ROC AUC (TPR vs FPR) (Alpha, Delta (1,2,3), BA.1)"
       , cex.axis = 1.5
       , cex.lab=1.5
       , main = "All values" 
       , cex.main = 1.5
       )
plot(  x = ROC_stats_filtered$w_all_ex_2_8_normMCC 
       , y = ROC_stats_filtered$w_all_ex_2_8_ROC_AUC
       , xlab = "Normalised MCC (Alpha, Delta (1,2,3), BA.1)"
       , ylab = "ROC AUC (TPR vs FPR) (Alpha, Delta (1,2,3), BA.1)"
       , cex.axis = 1.5
       , cex.lab=1.5
       , main = "ROC AUC > 0.7"
       , cex.main = 1.5
       #, main = "'True' positive EWS window: -30<= t <= +5" 
)

setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/non_TFPS/analysis/ROC_stats/" )
#' Choose filename
write.csv( ROC_stats_filtered , "ROC_stats_filtered_sd_non_tfps_5d.csv")
write.csv( ROC_stats_filtered , "ROC_stats_filtered_rt_non_tfps_5d.csv")
write.csv( ROC_stats_filtered , "ROC_stats_filtered_sd_non_tfps_10d.csv")
write.csv( ROC_stats_filtered , "ROC_stats_filtered_rt_non_tfps_10d.csv")

par (mfrow = c (3, 1))
# normMCC, ROC_AUC, F1, FMI
plot(  x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
     , y = ROC_stats_filtered$w_all_ex_2_8_normMCC
     , typ="l"
     , col = 1
     , main="ROC AUC > 0.7"
     , xlab = ""
     , ylab=""#normMCC"
     , ylim=c(0,1)
     ,cex.axis=1.5)
lines( x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
     , y = ROC_stats_filtered$w_all_ex_2_8_ROC_AUC
     , col=2 )
lines( x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
     , y = ROC_stats_filtered$w_all_ex_2_8_F1
     , col=3 )
lines( x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
     , y = ROC_stats_filtered$w_all_ex_2_8_FMI
     , col=4 )
legend("bottomright", legend=c("normMCC","ROC AUC","F1","FMI")
       , col=c(1,2,3,4)
       ,lty=c(1,1,1,1)
       ,cex=1.5)
# TP,FP,TN,FN
plot( x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
    , y = ROC_stats_filtered$w_all_ex_2_8_TP
    , col=5#"red"
    , typ="l"
    , xlab = ""
    , ylab="Number"
    , ylim=c(0,300)
    ,cex.axis=1.5)
lines( x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
     , y = ROC_stats_filtered$w_all_ex_2_8_FP
     , col=6#"orange"
     , typ="l" )
lines(  x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
        , y = ROC_stats_filtered$w_all_ex_2_8_TN 
        , col=7)
lines(  x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
        , y = ROC_stats_filtered$w_all_ex_2_8_FN
        , col=8)
legend("bottomright", legend=c("TP","FP","TN","FN")
       , col=c(5,6,7,8)
       ,lty=c(1,1,1,1)
       ,cex=1.5)
# TPR,FPR,TNR,PPV,NPV
plot(  x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
       , y = ROC_stats_filtered$w_all_ex_2_8_TPR
       , typ="l"
       ,col=9
       ,xlab = ""
       ,ylab=""#normMCC"
       ,ylim=c(0,1)
       ,cex.axis=1.5
)
lines(  x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
        , y = ROC_stats_filtered$w_all_ex_2_8_FPR
        , col=10)
lines(  x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
        , y = ROC_stats_filtered$w_all_ex_2_8_TNR
        , col=11)
lines(  x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
        , y = ROC_stats_filtered$w_all_ex_2_8_PPV
        , col=12)
lines(  x = seq( 1 , length( ROC_stats_filtered$w_all_ex_2_8_normMCC ) , 1 )
        , y = ROC_stats_filtered$w_all_ex_2_8_NPV
        , col=13)
legend("bottomright", legend=c("TPR aka sensitivity","FPR","TNR aka specificty","PPV aka precision","NPV (negative predictive value)")
       , col=c(9,10,11,12,13)
       ,lty=c(1,1,1,1,1)
       ,cex=1.5)

#' Search for leading indicators with a normMCC value for waves 3 to 7 
#' (because wave 2 impacted by not enough time before and wave 8 not a complete period after)
output_rt_normMCC = output_rt[ , c("w3_normMCC","w4_normMCC","w5_normMCC","w6_normMCC","w7_normMCC") ]
min_normMCC = data.frame( array( dim = c( dim( output_rt_normMCC )[ 1 ] , 1 ) ) )
for(i in 1 : nrow( output_rt_normMCC ) ){
  min_normMCC[ i , 1 ] = min( output_rt_normMCC[ i ,1:3 ] )
}
max(min_normMCC,na.rm=TRUE)
max_index = which( min_normMCC == max(min_normMCC,na.rm=TRUE) )
min_normMCC[max_index,]
View(t(output_rt[max_index,]))


#' Search for min value across various waves and then identify the leading indicator with the max of this min
output_rt_normMCC = output_rt[ , c("w3_normMCC","w4_normMCC","w7_normMCC") ]
min_normMCC = data.frame( array( dim = c( dim( output_rt_normMCC )[ 1 ] , 1 ) ) )
for(i in 1 : nrow( output_rt_normMCC ) ){
  min_normMCC[ i , 1 ] = min( output_rt_normMCC[ i ,1:3 ] )
}
max(min_normMCC,na.rm=TRUE)
max_index = which( min_normMCC == max(min_normMCC,na.rm=TRUE) )
min_normMCC[max_index,]
View(t(output_rt[max_index,]))



#' Checking for individual leading indicator set
li_index = which( names(li_df) == "CoMix_All_mean" ) -1 #' Subtract 1 as the first column is the date
ews_th_index = which( seq(-5,+5,0.05) == -0.75 ) 
seq(-5,+5,0.05)[ews_th_index]
conf = table( li_z_ews_rt_classification_array[ date_index_list[[ 9 ]] , li_index , ews_th_index ] )
TP = ifelse( is.na( conf["TP"] ) , 0 , conf["TP"] ) #' TP
FP = ifelse( is.na( conf["FP"] ) , 0 , conf["FP"] ) #' FP
TN = ifelse( is.na( conf["TN"] ) , 0 , conf["TN"] ) #' TN
FN = ifelse( is.na( conf["FN"] ) , 0 , conf["FN"] ) #' FN
TPR = TP / (TP + FN) #' TPR True Positive Rate aka Sensitivity, true positive rate, recall, hit rate
FPR = FP / (FP + TN) #' FPR False Positive Rate = 1 - Specificity
TNR = TN / (TN + FP) #' TNR = True Negative Rate aka Specificity, selectivity, true negative rate
PPV = TP / (TP + FP)  #' PPV aka Precision, positive predictive value
NPV = TN / (TN + FN) #' NPV = Negative predictive value
F1 = (2*TP) / (2*TP + FP + FN) #' F1
FMI = sqrt( PPV * TPR )  #' Fowlkes-Mallows index
#' Need to break down MCC calculation as otherwise results in NA with message 'NAs produced by integer overflow'
MCC = ((TP*TN) - (FP*FN)) / sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) #' cpmbines the 4 basic rates of confusion matrix
normMCC = (MCC + 1) /2
TP;FP;TN;FN;TPR;FPR;TNR;PPV;NPV;F1;FMI;MCC;normMCC

###########
#**OLD**
########################################
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