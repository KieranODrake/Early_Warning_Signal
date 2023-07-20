#' Create table of EWS results
#' Format: rows represent different leading indicators, columns represent the 7 
#' different EWS statistics for each wave, with one column for each wave showing
#'  the earliest EWS
#' 
#' #' create data frame to fill
#' ** only used first time created**
ews_results_df = data.frame( "leading_indicators" = "test"
                             #'wave 2
                         , "w2_start" = "test"
                         , "w2_EWS_earliest" = "test"
                         , "w2_EWS_diff" = "test"
                         , "w2_EWS_sd" = "test"
                         , "w2_EWS_skew" = "test"
                         , "w2_EWS_acf" = "test"
                         , "w2_EWS_Comp_sd_skew_acf" = "test"
                         , "w2_EWS_Comp_sd_skew" = "test"
                         , "w2_EWS_Comp_skew_acf" = "test"
                         , "w2_EWS_Comp_sd_acf" = "test"
                         #' wave 3
                         , "w3_start" = "test"
                         , "w3_EWS_earliest" = "test"
                         , "w3_EWS_diff" = "test"
                         , "w3_EWS_sd" = "test"
                         , "w3_EWS_skew" = "test"
                         , "w3_EWS_acf" = "test"
                         , "w3_EWS_Comp_sd_skew_acf" = "test"
                         , "w3_EWS_Comp_sd_skew" = "test"
                         , "w3_EWS_Comp_skew_acf" = "test"
                         , "w3_EWS_Comp_sd_acf" = "test"
                         #'wave 4
                         , "w4_start" = "test"
                         , "w4_EWS_earliest" = "test"
                         , "w4_EWS_diff" = "test"
                         , "w4_EWS_sd" = "test"
                         , "w4_EWS_skew" = "test"
                         , "w4_EWS_acf" = "test"
                         , "w4_EWS_Comp_sd_skew_acf" = "test"
                         , "w4_EWS_Comp_sd_skew" = "test"
                         , "w4_EWS_Comp_skew_acf" = "test"
                         , "w4_EWS_Comp_sd_acf" = "test"
                         #'wave 5
                         , "w5_start" = "test"
                         , "w5_EWS_earliest" = "test"
                         , "w5_EWS_diff" = "test"
                         , "w5_EWS_sd" = "test"
                         , "w5_EWS_skew" = "test"
                         , "w5_EWS_acf" = "test"
                         , "w5_EWS_Comp_sd_skew_acf" = "test"
                         , "w5_EWS_Comp_sd_skew" = "test"
                         , "w5_EWS_Comp_skew_acf" = "test"
                         , "w5_EWS_Comp_sd_acf" = "test"
                         #' wave 6
                         , "w6_start" = "test"
                         , "w6_EWS_earliest" = "test"
                         , "w6_EWS_diff" = "test"
                         , "w6_EWS_sd" = "test"
                         , "w6_EWS_skew" = "test"
                         , "w6_EWS_acf" = "test"
                         , "w6_EWS_Comp_sd_skew_acf" = "test"
                         , "w6_EWS_Comp_sd_skew" = "test"
                         , "w6_EWS_Comp_skew_acf" = "test"
                         , "w6_EWS_Comp_sd_acf" = "test"
                         #' wave 7
                         , "w7_start" = "test"
                         , "w7_EWS_earliest" = "test"
                         , "w7_EWS_diff" = "test"
                         , "w7_EWS_sd" = "test"
                         , "w7_EWS_skew" = "test"
                         , "w7_EWS_acf" = "test"
                         , "w7_EWS_Comp_sd_skew_acf" = "test"
                         , "w7_EWS_Comp_sd_skew" = "test"
                         , "w7_EWS_Comp_skew_acf" = "test"
                         , "w7_EWS_Comp_sd_acf" = "test"
                         #' wave 8
                         , "w8_start" = "test"
                         , "w8_EWS_earliest" = "test"
                         , "w8_EWS_diff" = "test"
                         , "w8_EWS_sd" = "test"
                         , "w8_EWS_skew" = "test"
                         , "w8_EWS_acf" = "test"
                         , "w8_EWS_Comp_sd_skew_acf" = "test"
                         , "w8_EWS_Comp_sd_skew" = "test"
                         , "w8_EWS_Comp_skew_acf" = "test"
                         , "w8_EWS_Comp_sd_acf" = "test"
)

#' Wave start dates (average where a wave has more than one start date as 
#' produced from the 12 models incorporating 'low resolution wave filter' - see
#' 21 June 2022 update slides)
wave_start_dates <- c( as.Date("")           #' 1 Wuhan
                       ,as.Date("2020-08-19") #' 2 B.1.177
                       ,as.Date("2020-11-29") #' 3 Alpha
                       ,as.Date("2021-05-11") #' 4 Delta
                       ,as.Date("2021-08-03") #' 5 Delta
                       ,as.Date("2021-09-27") #' 6 Delta
                       ,as.Date("2021-11-26") #' 7 Omicron
                       ,as.Date("2022-02-21") #' 8 Omicron
                       ,as.Date("")           #' 9 Omicron
)


#' Need to load ews_dates_array using 'EWS_plot.R'
#' loop through array and copy data to data frame
#' #' i is the
#' ** only use the next line the first time to set up the format**
#new_row = ews_results_df 
#'** after the first time use the following line to reset to NA**
new_row[,] = NA

for ( lead_ind in 1 : length( ews_names ) ) {
  new_row[, 1 ] = ews_names[ lead_ind ] #'"leading_indicators" = 
  for ( wave in 2 : 8 ){
    new_row[, (2 + ((wave-2)*10)) ] = wave_start_dates[ wave ]
    for ( ews_stat in 1 : 7 ){
      new_row[, (5 + ((wave-2)*10) + (ews_stat-1)) ] = ews_dates_array[ lead_ind , wave , ews_stat ]
      #new_row[, (4 + ((wave-2)*9) + (ews_stat-1)) ] = as.Date( ews_dates_array[ lead_ind , wave , ews_stat ] , origin="1970-01-01" ) 
    }
    #new_row[ , (3 + ((wave-2)*9)) ] = min( new_row[ (4+((wave-2)*9)) : (10+((wave-2)*9)) ] , na.rm = FALSE )
    temp_stats =  new_row[ (5+((wave-2)*10)) : (11+((wave-2)*10)) ]
    new_row[ , (3 + ((wave-2)*10)) ] = ifelse( length( temp_stats[ !is.na( temp_stats ) ] ) == 0 , NA_real_ , min( temp_stats , na.rm = TRUE ))
    for ( item in temp_stats ){
      item = as.Date(item , origin = "1970-01-01")
    }
    new_row[ , (3 + ((wave-2)*10)) ] = as.Date( new_row[ , (3 + ((wave-2)*10)) ] , origin = "1970-01-01" )
    new_row[ , (4 + ((wave-2)*10)) ] = new_row[ , (3 + ((wave-2)*10)) ] - new_row[ , (2 + ((wave-2)*10)) ]
  }
  ews_results_df = rbind( ews_results_df , new_row )
}

#' Remove first row which was the test in creating the data frame
#' #**only necessary for the time**
ews_results_df = ews_results_df[-1 , ] 


View(ews_results_df)

#' Just the EWS lead/lag. Note that -ve is an EWS before the wave start (inflection) and +ve is after
ews_diff_df = ews_results_df[,c(1,4,14,24,34,44,54,64)]
View(ews_diff_df)

#' Convert numbers to dates
for ( item in ews_results_df[2:nrow(ews_results_df),2:71] ){ 
  #if (!is.na( item )){ item = as.Date( as.numeric( item ) , origin = "1970-01-01") }
  item = as.Date( as.numeric( item ) , origin = "1970-01-01")
}
 
#' Save data frames
setwd('C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Analysis')
saveRDS( ews_results_df , file = "EWS_dates_2022_10_14.rds" )
saveRDS( ews_diff_df , file = "EWS_lead_lag_2022_10_14.rds" )

a = readRDS("EWS_dates.rds" )
View(a)
ews_diff_df = readRDS("EWS_lead_lag_2022_10_14.rds")
View(ews_diff_df)
ews_results_df = readRDS("EWS_dates_2022_10_14.rds")
View(ews_results_df)
#' Relable rows from 1 to n
library(data.table)
(setattr(ews_results_df, "row.names", c(as.character(seq(1,496,1)))))

####################################
#' plot data
#'

#' Change names of Google Mobility as too long
ews_diff_df[58,1] = "Google_mobility_retail_recreation"
ews_diff_df[59,1] = "Google_mobility_grocery_pharmacy"
ews_diff_df[60,1] = "Google_mobility_parks"
ews_diff_df[61,1] = "Google_mobility_transit_stations"
ews_diff_df[62,1] = "Google_mobility_workplaces"
ews_diff_df[63,1] = "Google_mobility_residential"
ews_diff_df[64,1] = "Google_mobility_inv_parks"

#' what is the min and max lead/lag?
min( as.numeric( ews_diff_df[ , 2:8 ] ) , na.rm = TRUE )

#' plot other leading indicators
par(mar=c(17,6,5,2))
par(mfrow=c(1,1))
par(yaxt="s")

row_min = 1 #1 #79
row_max = 64 #78 #140
plot(ews_diff_df[row_min:row_max,2], axes=FALSE, xlab="",ylab="EWS lead (-ve) / lag (+ve)",ylim=c(-100,100),pch=16,xlim = c(0,row_max+1))
#plot((ews_diff_df[,2] - shift_df[,2]), axes=FALSE, xlab="",ylab="EWS lead (-ve) / lag (+ve)",ylim=c(-100,100),pch=16)
#axis(2)
axis(1, at=seq_along(ews_diff_df[row_min:row_max,1]),labels=as.character(ews_diff_df[row_min:row_max,1]), las=2, cex.axis=0.75)
axis(side=2, yaxp=c(-100,100,20), las=2, cex.axis=0.75,tick=TRUE,lty="solid",col="black")
#' Add grid to readability
abline(v=c(seq(row_min,row_max,1)),col="gray")
abline(h=c(seq(-100,100,10)),col="gray")
abline(h=0)
box()
points(ews_diff_df[row_min:row_max,3],col="red",pch=15)
points(ews_diff_df[row_min:row_max,4],col="blue",pch=17)
points(ews_diff_df[row_min:row_max,5],col="green",pch=18)
points(ews_diff_df[row_min:row_max,6],col="darkgreen",pch=19)
points(ews_diff_df[row_min:row_max,7],col="orange",pch=20)
points(ews_diff_df[row_min:row_max,8],col="purple",pch=1)
legend("topright", legend=c("B.1.177","Alpha","Delta 1","Delta 2","Delta 3","Omicron 1", "Omicron 2")
       , pch=c(16,15,17,18,19,20,1)
       , col=c("black","red","blue","green","darkgreen","orange","purple"))

#' plot TFP Scan leading indicators
par(mar=c(17,6,5,2))
par(mfrow=c(1,1))
par(yaxt="s")

row_min = 209 #65 #137 #209 #233 #305 #377 #449
row_max = 232 #136 #208 #232 #304 #376 #448 449+3*24-1
plot(ews_diff_df[row_min:row_max,2], axes=FALSE, xlab="",ylab="EWS lead (-ve) / lag (+ve)",ylim=c(-100,100),pch=16)#xlim = c(0,row_max+1))
#plot((ews_diff_df[,2] - shift_df[,2]), axes=FALSE, xlab="",ylab="EWS lead (-ve) / lag (+ve)",ylim=c(-100,100),pch=16)
#axis(2)
axis(1, at=seq_along(ews_diff_df[row_min:row_max,1]),labels=as.character(ews_diff_df[row_min:row_max,1]), las=2, cex.axis=0.75)
axis(side=2, yaxp=c(-100,100,20), las=2, cex.axis=0.75,tick=TRUE,lty="solid",col="black")
#' Add grid to readability
abline(v=c(seq(row_min-(row_min-1),row_max-(row_min-1),1)),col="gray")
abline(h=c(seq(-100,100,10)),col="gray")
abline(h=0)
box()
points(ews_diff_df[row_min:row_max,3],col="red",pch=15)
points(ews_diff_df[row_min:row_max,4],col="blue",pch=17)
points(ews_diff_df[row_min:row_max,5],col="green",pch=18)
points(ews_diff_df[row_min:row_max,6],col="darkgreen",pch=19)
points(ews_diff_df[row_min:row_max,7],col="orange",pch=20)
points(ews_diff_df[row_min:row_max,8],col="purple",pch=1)
legend("topright", legend=c("B.1.177","Alpha","Delta 1","Delta 2","Delta 3","Omicron 1", "Omicron 2")
       , pch=c(16,15,17,18,19,20,1)
       , col=c("black","red","blue","green","darkgreen","orange","purple"))

#'set number of days to shift the number of days between EWS and wave start by 
shift_days = 14
shift_df = data.frame( ews_diff_df[ , 1 ]
                       , as.numeric( ews_diff_df[ , 2 ] ) - shift_days
                       , as.numeric( ews_diff_df[ , 3 ] ) - shift_days
                       , as.numeric( ews_diff_df[ , 4 ] ) - shift_days
                       , as.numeric( ews_diff_df[ , 5 ] ) - shift_days
                       , as.numeric( ews_diff_df[ , 6 ] ) - shift_days
                       , as.numeric( ews_diff_df[ , 7 ] ) - shift_days
                       , as.numeric( ews_diff_df[ , 8 ] ) - shift_days
)

#' Plot just for B.1.177 with shift 
plot_col = "purple" #"black"
pch_val = 21 #16
plot_column = 8 #2
#' vlgr etc
row_min = 1 #1 #79
row_max = 78 #78 #140
plot(ews_diff_df[row_min:row_max,plot_column], axes=FALSE, xlab="",ylab="EWS lead (-ve) / lag (+ve)",ylim=c(-100,100),xlim = c(0,row_max+1),col=plot_col,pch=pch_val,cex=2)
points(shift_df[row_min:row_max,plot_column],pch=pch_val,col="grey",cex=2)
#' Add grid to readability
abline(v=c(seq(row_min,row_max,1)),col="gray")
abline(h=c(seq(-100,100,10)),col="gray")
abline(h=0)
axis(1, at=seq_along(ews_diff_df[row_min:row_max,1]),labels=as.character(ews_diff_df[row_min:row_max,1]), las=2, cex.axis=0.75)
axis(side=2, yaxp=c(-100,100,20), las=2, cex.axis=0.75,tick=TRUE,lty="solid",col="black")
box()
#' Plot just for B.1.177 with shift : other leading indicators
row_min = 79 #1 #79
row_max = 142 #78 #142
plot(ews_diff_df[row_min:row_max,plot_column], axes=FALSE, xlab="",ylab="EWS lead (-ve) / lag (+ve)",ylim=c(-100,100),pch=pch_val,col=plot_col,cex=2)
points(shift_df[row_min:row_max,plot_column],pch=pch_val,col="grey",cex=2)
#' Add grid to readability
abline(v=c(seq(row_min-78,row_max-78,1)),col="gray")
abline(h=c(seq(-100,100,10)),col="gray")
abline(h=0)
axis(1, at=seq_along(ews_diff_df[row_min:row_max,1]),labels=as.character(ews_diff_df[row_min:row_max,1]), las=2, cex.axis=0.75)
axis(side=2, yaxp=c(-100,100,20), las=2, cex.axis=0.75,tick=TRUE,lty="solid",col="black")
box()

points(ews_diff_df[row_min:row_max,3],col="red",pch=15)
points(ews_diff_df[row_min:row_max,4],col="blue",pch=17)
points(ews_diff_df[row_min:row_max,5],col="green",pch=18)
points(ews_diff_df[row_min:row_max,6],col="darkgreen",pch=19)
points(ews_diff_df[row_min:row_max,7],col="orange",pch=20)
points(ews_diff_df[row_min:row_max,8],col="purple",pch=1)
legend("topright", legend=c("B.1.177","Alpha","Delta 1","Delta 2","Delta 3","Omicron 1", "Omicron 2")
       , pch=c(16,15,17,18,19,20,1)
       , col=c("black","red","blue","green","darkgreen","orange","purple"))

######################
#' Calculate the best EWS for each set of leading indicators
li_type = c("test_hosp_shift","Ct_value","positivity_rate","CoMix","Google_mobility","vlgr","vlgr/mean_cluster_size","lgr_max","lgr_mean","lgr_wtd_mean","clock_outlier_max","clock_outlier_mean","clock_outlier_z_x_lgr_z")
library(data.table)
(setattr(ews_diff_df, "row.names", c(as.character(seq(1,520,1)))))
View(ews_diff_df)
row_s = c(1, 5,43,45,58,65 ,209,233,305,377,449,473,497)
row_e = c(4,42,44,57,64,208,232,304,376,448,472,496,520)
ews_diff_earliest_df = ews_diff_df[1:length(li_type),]
ews_diff_earliest_df[ , ] = NA
ews_diff_earliest_df[ 1:length(li_type) , 1 ] = li_type

#' Cycle through waves
for (i in 1:length( li_type )){
  for (j in 2:8 ){
    ews_diff_earliest_df[ i , j ] = ifelse( all( is.na( ews_diff_df[ row_s[i]:row_e[i] , j ] ) )
                                            , NA
                                            , min( as.numeric( ews_diff_df[ row_s[i]:row_e[i] , j ] ) , na.rm = TRUE )
                                            )
  }
}

View(ews_diff_earliest_df)

par(mar=c(17,6,5,2))
par(mfrow=c(1,1))
par(yaxt="s")

row_min = 1 
row_max = 12 
plot(ews_diff_earliest_df[row_min:row_max,2], axes=FALSE, xlab="",ylab="EWS lead (-ve) / lag (+ve)",ylim=c(-100,100),pch=16,xlim = c(0,row_max+1),cex=2)
#plot((ews_diff_df[,2] - shift_df[,2]), axes=FALSE, xlab="",ylab="EWS lead (-ve) / lag (+ve)",ylim=c(-100,100),pch=16)
#axis(2)
axis(1, at=seq_along(ews_diff_earliest_df[row_min:row_max,1]),labels=as.character(ews_diff_earliest_df[row_min:row_max,1]), las=2, cex.axis=1)
axis(side=2, yaxp=c(-100,100,20), las=2, cex.axis=0.75,tick=TRUE,lty="solid",col="black")
#' Add grid to readability
abline(v=c(seq(row_min,row_max,1)),col="gray")
abline(h=c(seq(-100,100,10)),col="gray")
abline(h=0)
box()
points(ews_diff_earliest_df[row_min:row_max,3],col="red",pch=15,cex=2)
points(ews_diff_earliest_df[row_min:row_max,4],col="blue",pch=17,cex=2)
points(ews_diff_earliest_df[row_min:row_max,5],col="green",pch=18,cex=2)
points(ews_diff_earliest_df[row_min:row_max,6],col="darkgreen",pch=19,cex=2)
points(ews_diff_earliest_df[row_min:row_max,7],col="orange",pch=20,cex=2)
points(ews_diff_earliest_df[row_min:row_max,8],col="purple",pch=1,cex=2)
legend("topright", legend=c("B.1.177","Alpha","Delta 1","Delta 2","Delta 3","Omicron 1", "Omicron 2")
       , pch=c(16,15,17,18,19,20,1)
       , col=c("black","red","blue","green","darkgreen","orange","purple"))
