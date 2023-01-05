#' Calculate correlation between variance of cluster logistic growth rates 
#' from TFP scanner and Rt 
#' Author: Kieran Drake
#' Date: December 2022

library(data.table)
library(gratia)
library(mgcv)
library(lubridate)
library(sqldf)
library(zoo)

####    Load Rt data    ####
#' Rt data downloaded from https://www.gov.uk/guidance/the-r-value-and-growth-rate 221207_R_and_growth_rate_time_series_for_publication_v1.0.ods
#' Asterisks removed from some dates in date column and column reformatted to 'Date' in Excel
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_12/analysis/correlation") 
rt = read.csv("Rt.csv")
rt$mid <- ( rt$lower + rt$upper ) / 2
rt$date <- as.Date( rt$date , format = "%d/%m/%Y") # Change format of date
rt$time <- lubridate::decimal_date( rt$date ) # Add column for decimal date

####    Calculate derivative of Rt (dRt/dt)    ####
deg_free_k = 100
GAM_smooth_function = "tp"
m <- mgcv::gam( mid ~ s( time, bs=GAM_smooth_function, k = deg_free_k), data = rt )
significance_level = 0.95
drt_dt = gratia::derivatives( m , term="s(time)" , n = length( rt$mid ) , level = significance_level )#, partial_match = TRUE )
#'plot
par(mfrow=c(2,1))
plot(rt$time,rt$mid,typ="l",lwd=2,col="blue",xlab="year",ylab="Rt",ylim=c(0.5,2.5))
lines( rt$time,rt$lower , typ="l" , col = "blue" , lwd=1, lty=2 )
lines( rt$time,rt$upper , typ="l" , col = "blue" , lwd=1, lty=2 )
lines( m$model[[2]] , m$fitted.values , typ="l" , col = "red" , lwd=2 )
legend("topright",legend=c("Rt (midpoint of 90% CI for England)", "Rt 90% CI band",paste0("Rt model fit (spline=",GAM_smooth_function,", k=",deg_free_k,")"),"Rt = 1")
        , lty=c(1,2,1,2), lwd=c(2,1,2,2),col=c("blue","blue","red","black")
       )
abline(h=1,lty=2)

plot(drt_dt$data,drt_dt$derivative,typ="l",xlab="time",ylab="d(Rt)/dt")
abline(h=0)

####    Load Var(LGR) data    ####
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_12/analysis/correlation") 
#' There are six different types of Var(LGR) statistics from TFP scanner...
var_types = c(  "tfps_vlgr_gam_pop.csv"
              , "tfps_vlgr_gam_samp.csv"
              , "tfps_vlgr_pop.csv"
              , "tfps_vlgr_samp.csv"
              , "tfps_vlgr_simple_pop.csv"
              , "tfps_vlgr_simple_samp.csv"
              )
#' ... three types of p-value filter...
pval = c("no","001","005")
#' ... and 10 different LGR threshold levels for replacing sub-clusters with parent clusters
lgr_th = c( "false" , "060", "065", "070", "075", "080", "085", "090", "095", "100" )
#' Calculate expected number of parameter combinations
par_exp = length(var_types) * length(pval) * length(lgr_th) * 24 #' 24 is the number of columns in each file and the number of combinations of min cluster age, max cluster age and min descendants
#' Initialise variables
file_list = list() ; counter = 0 ; counter2 = 0
#' Loop through all files and compile data into a single dataframe
for ( pv in pval ){
 for ( lgr in lgr_th ){
  for( vt in var_types ){
    counter = counter + 1
    file_list[counter] = paste0("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_11/analysis_new_pango_stat_2/large_cluster_adjust_lgr_threshold_",lgr,"/p_val_filter_" , pv , "/dataframes_statistics/" , vt )
    #' Copy all columns (representing different combinations of min cluster age, max cluster age and min descendants) from the file to the new dataframe
    for( i in 2:25 ){
      counter2 = counter2 + 1
      temp_lead_ind_stats = data.table::fread( paste0(file_list[ counter ]) ) ; temp_lead_ind_stats = as.data.frame(temp_lead_ind_stats) ; temp_lead_ind_stats[,1] = NULL #' Need to remove the first column as it just contains row numbes in the .csv files
      #' Create dataframe and add dates from the first file (which should be the same for all subsequent files)
      if( counter2 == 1 ){
        vlgr_data_df = data.frame( temp_lead_ind_stats[ , 1 ] )  
      }
      vlgr_data_df[ , i + (24*(counter-1)) ] = temp_lead_ind_stats[ , i ]
      colnames( vlgr_data_df )[ i + (24*(counter-1)) ] <- paste0(    strsplit( vt , split = ".csv" ) , "_"
                                                                    , lgr , "_"
                                                                    ,  pv , "_"
                                                                    , colnames( temp_lead_ind_stats )[ i ]
                                                                 )
      message(counter2, " of ", par_exp,". ", round(100*counter2/par_exp,1),"%. Just completed: p-value ",pv,", parent/sub LGR threshold ",lgr,", leading indicator ", strsplit( vt , split = ".csv" ), " and scan parameters ", colnames( temp_lead_ind_stats )[ i ] )
    }
    #message(counter2, " of ", par_exp,". ", round(100*counter2/par_exp,1),"%. Just completed: p-value ",pv,", parent/sub LGR threshold ",lgr,", and leading indicator ", strsplit( vt , split = ".csv" ) )
  }
 }
}

####    Calculate correlation, including with time shifts   ####
#' Merge d(Rt)/dt data with an empty dataframe with the leading indicator dates.
#' Then use xts::na.approx to fill in the empty dates
#' Then pull out the values for the leading indicator dates
#' Then calculate correlations between d(Rt)/dt and leading indicator values

#' create dates required to calculate correlation 
dates_extended <- data.frame( "time" = lubridate::decimal_date( vlgr_data_df[,1] ) , "derivative" = NA)
drt_dt_trimmed <- as.data.frame( drt_dt[,c("data","derivative")] )
names( drt_dt_trimmed )[c(1,2)] <- c("time","derivative")
dates_extended <- rbind( dates_extended , drt_dt_trimmed )
dates_extended <- data.table::setorder( dates_extended , cols = "time" )
#' Fill in d(Rt)/dt for dates with no values
dates_extended[,3] <- na.approx(dates_extended[,2])
dates_extended[,4] <- na.approx(dates_extended[,2] , 1:400 )
dates_extended[,5] <- na.spline(dates_extended[,2] )
dates_extended[,6] <- na.spline(dates_extended[,2] , 1:400 )
dates_extended[,7] <- rowMeans(dates_extended[,c(3,4,5,6)])
#' Final dataframe of d(Rt)/dt to compare against leading indicator data  
lead_ind_dates <- data.frame( "date" = vlgr_data_df[,1] , "time" = lubridate::decimal_date( vlgr_data_df[,1] ) , NA )
drt_dt_extended <- subset( dates_extended , dates_extended$time %in% lead_ind_dates$time )
#' Checking interpolation for d(Rt)/dt
par(mfrow = c(2,3))
plot(drt_dt$data,drt_dt$derivative,typ="l",xlab="time",ylab="d(Rt)/dt",cex.axis=1.5,cex.lab=1.5)#,col="blue")
lines(drt_dt_extended$time,drt_dt_extended[,3],col="red"); abline(h=0)
plot(drt_dt$data,drt_dt$derivative,typ="l",xlab="time",ylab="d(Rt)/dt",cex.axis=1.5,cex.lab=1.5)#,col="blue")
lines(drt_dt_extended$time,drt_dt_extended[,4],col="green"); abline(h=0)
plot(drt_dt$data,drt_dt$derivative,typ="l",xlab="time",ylab="d(Rt)/dt",cex.axis=1.5,cex.lab=1.5)#,col="blue")
lines(drt_dt_extended$time,drt_dt_extended[,5],col="blue"); abline(h=0)
plot(drt_dt$data,drt_dt$derivative,typ="l",xlab="time",ylab="d(Rt)/dt",cex.axis=1.5,cex.lab=1.5)#,col="blue")
lines(drt_dt_extended$time,drt_dt_extended[,6],col="orange"); abline(h=0)
plot(drt_dt$data,drt_dt$derivative,typ="l",xlab="time",ylab="d(Rt)/dt",cex.axis=1.5,cex.lab=1.5)#,col="blue")
lines(drt_dt_extended$time,drt_dt_extended[,7],col="violet"); abline(h=0)
#' the mean of the four different na.approx() seems to be the best fit to the derivative: drt_dt_extended[,7]
par(mfrow = c(1,1))
plot(drt_dt$data,drt_dt$derivative^3,typ="l",xlab="time",ylab="d(Rt)/dt",cex.axis=1.5,cex.lab=1.5,lwd=2)
lines(drt_dt_extended$time,drt_dt_extended[,3]^3,col="red",lwd=2)
lines(drt_dt_extended$time,drt_dt_extended[,4]^3,col="green",lwd=2)
lines(drt_dt_extended$time,drt_dt_extended[,5]^3,col="blue",lwd=2)
lines(drt_dt_extended$time,drt_dt_extended[,6]^3,col="orange",lwd=2)
lines(drt_dt_extended$time,drt_dt_extended[,7]^3,col="violet",lwd=2)

par(mfrow = c(1,4))
plot(drt_dt$data,drt_dt$derivative^3,typ="l",xlab="time",ylab="(d(Rt)/dt)^3",cex.axis=1.5,cex.lab=1.5)#,col="blue")
lines(drt_dt_extended$time,drt_dt_extended[,3]^3,col="red"); abline(h=0)
plot(drt_dt$data,drt_dt$derivative^3,typ="l",xlab="time",ylab="(d(Rt)/dt)^3",cex.axis=1.5,cex.lab=1.5)#,col="blue")
lines(drt_dt_extended$time,drt_dt_extended[,4]^3,col="green"); abline(h=0)
plot(drt_dt$data,drt_dt$derivative^3,typ="l",xlab="time",ylab="(d(Rt)/dt)^3",cex.axis=1.5,cex.lab=1.5)#,col="blue")
lines(drt_dt_extended$time,drt_dt_extended[,5]^3,col="blue"); abline(h=0)
plot(drt_dt$data,drt_dt$derivative^3,typ="l",xlab="time",ylab="(d(Rt)/dt)^3",cex.axis=1.5,cex.lab=1.5)#,col="blue")
lines(drt_dt_extended$time,drt_dt_extended[,6]^3,col="orange"); abline(h=0)
par(mfrow = c(1,1))
plot(drt_dt$data,drt_dt$derivative^3,typ="l",xlab="time",ylab="(d(Rt)/dt)^3",cex.axis=1.5,cex.lab=1.5)#,col="blue")
lines(drt_dt_extended$time,drt_dt_extended[,7]^3,col="violet"); abline(h=0)



#' This doesn't seem to work
a = as.vector(PredictMat( m$smooth[[1]], data = lead_ind_dates) %*% coef(m)[ grepl( names(coef(m)), patt = 's\\(time\\)' ) ] )
plot(drt_dt$data,drt_dt$derivative,typ="l",xlab="time",ylab="d(Rt)/dt")
lines(lead_ind_dates$time,a,typ="l",col="red")

#' Test plot to show leading indicator and d(Rt)/dt
par(mfrow = c(1,1))
plot(lubridate::decimal_date(vlgr_data_df[,1]),vlgr_data_df[,2],ylim=c(0,10),typ="l")
lines(drt_dt_extended$time,drt_dt_extended[,7]*0.1+1,col="red",typ="l")

#' Compute correlations across all Var(LGR) leading indicators
#' Also shift up to 28 days (moving leading indicator date forward in time)
shift <- function( x , n ){
  c( x[ -( seq( n ) ) ], rep( NA , n ) )
}

#' alternative to try to shift both ways
#shift <- function( x , n ){
#  ifelse( (n = 0) | (n > 0)
#          , c( x[ -( seq( n ) ) ], rep( NA , n ) )
#          , c( rep( NA , -n ), x[ +( seq( -n ) ) ], rep( NA , n ) )
#          )
#}

#' shift leading indicator data backwards in time #** Not sure this is the correct shift direction - higher Var(LGR) means new higher transmissible variant and so Rt with increase so higher d(Rt)/dt**
cor_li_drt_dt_df = data.frame()
for (k in 0:28){ #seq(-28,28,1)){
  for(j in 2:ncol(vlgr_data_df) ){
    if( k == 0 ){ #' i.e. no shift because shift() doesn't work as might be expected for n=0
      vlgr_shift = vlgr_data_df[ , j ] 
    } else { 
      vlgr_shift = shift( vlgr_data_df[ , j ] , k )
    }
    #vlgr_shift = shift( vlgr_data_df[,j] , k )
    #pre_cor_df = data.frame( vlgr_data_df[,j] , drt_dt_extended[,7] )
    pre_cor_df = data.frame( vlgr_shift , drt_dt_extended[ 1:length(vlgr_shift) , 7 ] )
    pre_cor_df = subset( pre_cor_df , !is.na(pre_cor_df[,1]) & !is.na(pre_cor_df[,2]) )
    cor_li_drt_dt_df[ k+1 , j-1 ] = cor( pre_cor_df[ , 1 ] , pre_cor_df[ , 2 ] )
    message("Shift ",k,", vlgr ",j)
  }
}

#' shift leading indicator data backwards in time 
cor_li_drt_dt_df_2 = data.frame()
for (k in 0:28){ #seq(-28,28,1)){
  for(j in 2:ncol(vlgr_data_df) ){
    if( k == 0 ){ #' i.e. no shift because shift() doesn't work as might be expected for n=0
      drt_dt_extended_shift = drt_dt_extended[ , 7 ]
    } else{
      drt_dt_extended_shift = shift( drt_dt_extended[ , 7 ] , k )
    }
    #vlgr_shift = shift( vlgr_data_df[,j] , k )
    #drt_dt_extended_shift = shift( drt_dt_extended[,7] , k )
    #pre_cor_df = data.frame( vlgr_shift , drt_dt_extended[ 1:length(vlgr_shift) , 7 ] )
    pre_cor_df = data.frame( vlgr_data_df[ 1:length(drt_dt_extended_shift) , j ] , drt_dt_extended_shift )
    pre_cor_df = subset( pre_cor_df , !is.na(pre_cor_df[,1]) & !is.na(pre_cor_df[,2]) )
    cor_li_drt_dt_df_2[ k+1 , j-1 ] = cor( pre_cor_df[ , 1 ] , pre_cor_df[ , 2 ] )
    message("Shift ",k,", vlgr ",j)
  }
}

#' check zero shift in both directions is equal
cor_li_drt_dt_df[1,1:10] == cor_li_drt_dt_df_2[1,1:10]

max(cor_li_drt_dt_df_2 , na.rm=TRUE)
max(cor_li_drt_dt_df_2[1,] , na.rm=TRUE)
max(cor_li_drt_dt_df , na.rm=TRUE)
max(cor_li_drt_dt_df[1,] , na.rm=TRUE)

plot(cor_li_drt_dt_df[,1])
plot(cor_li_drt_dt_df_2[,1])

#' Plot surface of time shifted correlations
library(plotly)
#cor_li_drt_dt_df_combined = rbind( rev( cor_li_drt_dt_df_2[2:29,1:2] ) , cor_li_drt_dt_df[,])
cor_li_drt_dt_matrix = as.matrix(cor_li_drt_dt_df,nrow=nrow(cor_li_drt_dt_df),ncol=ncol(cor_li_drt_dt_df))
#cor_li_drt_dt_matrix = as.matrix( cor_li_drt_dt_df_combined, nrow=nrow(cor_li_drt_dt_df_combined), ncol=ncol(cor_li_drt_dt_df_combined))
fig <- plot_ly(z = ~cor_li_drt_dt_matrix) #volcano)
fig <- fig %>% add_surface()
fig <- fig %>% layout(title = 'Manually Specified Labels'
                      ,  xaxis = list(title = 'leading indicator parameter set')
                      ,  yaxis = list(title = 'Number of days Var(LGR) shifted relative to d(Rt)/dt')
                      #,  zaxis = list(title = 'Pearson correlation coefficient')
                      , legend = list(title=list(text='<b> Pearson correlation coefficient </b>')))
fig
plot_ly(z = cor_li_drt_dt_matrix , type = "heatmap", x = seq(1,4320,1), y = seq(0,28,1) )

cor_li_drt_dt_matrix_2 = as.matrix(cor_li_drt_dt_df_2,nrow=nrow(cor_li_drt_dt_df_2),ncol=ncol(cor_li_drt_dt_df_2))
fig <- plot_ly(z = ~cor_li_drt_dt_matrix) #volcano)
fig <- fig %>% add_surface()
fig

#' Compute mean, median and maximum correlations by number of days shifted across all Var(LGR) parameter sets
for (i in 1:nrow(cor_li_drt_dt_df)){
  print( paste( i-1
                , mean(   as.numeric( cor_li_drt_dt_df[ i, ] ) , na.rm=TRUE )
                , median( as.numeric( cor_li_drt_dt_df[ i, ] ) , na.rm=TRUE )
                , max( as.numeric( cor_li_drt_dt_df[ i, ] ) , na.rm=TRUE )
               ) 
        )
}

#' Violin plots to show distribution of Pearson Correlation Coefficients across the nuber of days shifted
library(plotly)

df <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/violin_data.csv")

fig <- df %>%
  plot_ly(
    x = ~day,
    y = ~total_bill,
    split = ~day,
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
      title = "Day"
    ),
    yaxis = list(
      title = "Total Bill",
      zeroline = F
    )
  )

fig

#' Search for parameter set with highest correlation
#' Highest correlation with no shift
j = which(cor_li_drt_dt_df[1,] == max( as.numeric( cor_li_drt_dt_df[1,] ) , na.rm=TRUE ) )
cor_li_drt_dt_df[1,j]
names(vlgr_data_df)[j+1]

#' Highest correlation with any shift
for (i in 1:nrow(cor_li_drt_dt_df) ){ 
  max_1 = max( as.numeric( cor_li_drt_dt_df[i,] ) , na.rm=TRUE )
  t1 = which(cor_li_drt_dt_df[i,] == max_1 ) 
  max_2 = max( as.numeric( cor_li_drt_dt_df_2[i,] ) , na.rm=TRUE )
  t2 = which(cor_li_drt_dt_df_2[i,] == max_2 )
  print( paste( i , t1 , round( max_1 , 4 ) , t2 , round( max_2 , 4 ) ) )
}
  cor_li_drt_dt_df[1,j]
  names(vlgr_data_df)[j+1]


par(mfrow = c(1,1))
plot(lubridate::decimal_date(vlgr_data_df[,1]),vlgr_data_df[,j+1],ylim=c(0,10),typ="l",main=names(vlgr_data_df)[j+1],ylab="",xlab="year")
lines(drt_dt_extended$time,drt_dt_extended[,7],col="red",typ="l")
par(mfrow = c(2,1))
plot(drt_dt_extended$time,drt_dt_extended[,7],col="red",typ="l",main=names(vlgr_data_df)[j+1],ylab="",xlab="year",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(lubridate::decimal_date(vlgr_data_df[,1]),vlgr_data_df[,j+1]*10)
legend("bottomleft",legend=c("Var(LGR) x10","d(Rt)/dt"),col=c("black","red"),lty=c(1,1),cex=1.5)
cor(drt_dt_extended[,7] , vlgr_data_df[,j+1] ,na.rm=TRUE)

plot(lubridate::decimal_date(vlgr_data_df[,1]),vlgr_data_df[,j+1],typ="l",main=names(vlgr_data_df)[j+1],ylab="",xlab="year",ylim=c(0,2),cex.axis=1.5,cex.lab=1.5)
lines(rt$time,rt$mid,col="blue")
abline(h=1,lty=2,col="blue")
legend("topright",legend=c("Var(LGR)","Rt","Rt = 1"),col=c("black","blue","blue"),lty=c(1,1,2),cex=1.5)

par(mfrow=c(1,1))
lead_ind_index = j
plot(  x = seq(-28,28,1)
       , y = c( rev( cor_li_drt_dt_df_2[2:29,lead_ind_index] ) , cor_li_drt_dt_df[,lead_ind_index])
       , xlab = "Number of days d(Rt)/dt shifted relative to Var(LGR)", ylab="Pearson correlation coefficient"
       ,cex.axis=1.5, cex.lab=1.5)
lines(  x = seq(-28,28,1)
        , y = c( rev( cor_li_drt_dt_df_2[2:29,lead_ind_index] ) , cor_li_drt_dt_df[,lead_ind_index])
)
abline(v=0,h=0)

#' search for parameter sets with correlation above x
k = which(cor_li_drt_dt_df[1,] > 0.4 )
cor_li_drt_dt_df[1,k]
names(vlgr_data_df)[k+1]

hist(as.numeric(cor_li_drt_dt_df[29,]))

par(mfrow=c(3,10))
for (i in 1:29){
  hist(as.numeric(cor_li_drt_dt_df[i,]))
}
lapply( cor_li_drt_dt_df[ seq(1,29,1) , ]  , FUN=hist )


####    Calculate rolling correlation - different time periods     ####
#'**perhaps try 7, 14, 28 day rolling correlation, then make available for EV to look at ?**
#' Need to extrapolate to daily data so probably only look at minimum 28 day rolling correlation to reduce % that are synthetic data points
#' Create all calendar dates between start and end of Var(LGR) data

calendar_dates = data.frame( as.Date( seq( min(vlgr_data_df[,1]) , max(vlgr_data_df[,1]) , 1 ) ) ) 
calendar_dates[,2] = lubridate::decimal_date(calendar_dates[,1])
names(calendar_dates) <- c("date","time")
#' Select dates with no Var(LGR) data and set data points as NA
missing_dates = subset(calendar_dates[,1],  !(calendar_dates[,1] %in% vlgr_data_df[,1]) )
empty_mat = matrix( NA , nrow = length( missing_dates ) , ncol = ncol( vlgr_data_df ) - 1 )
missing_dates = cbind( missing_dates , as.data.frame( empty_mat ) )
names(missing_dates) <- names(vlgr_data_df)
#' Merge the dates with and without Var(LGR) data
vlgr_data_full_calendar = as.data.frame( rbind( vlgr_data_df , missing_dates) )
names(vlgr_data_full_calendar)[1] = "date"
vlgr_data_full_calendar <- data.table::setorder( vlgr_data_full_calendar , cols = "date" )

#' Interpolate Var(LGR) for dates with no values
vlgr_data_full_calendar_interpolate_1 = vlgr_data_full_calendar
vlgr_data_full_calendar_interpolate_2 = vlgr_data_full_calendar
vlgr_data_full_calendar_interpolate_3 = vlgr_data_full_calendar
vlgr_data_full_calendar_interpolate_4 = vlgr_data_full_calendar

for (m in 2:ncol(vlgr_data_full_calendar_interpolate_1)){
  #' linear interpolation
  vlgr_data_full_calendar_interpolate_1[,m] = zoo::na.approx( vlgr_data_full_calendar[,m], na.rm=FALSE )
  vlgr_data_full_calendar_interpolate_2[,m] = zoo::na.approx( vlgr_data_full_calendar[,m], na.rm=FALSE, 1:nrow(vlgr_data_full_calendar) )
  #' Cubic spline interpolation
  vlgr_data_full_calendar_interpolate_3[,m] = zoo::na.spline( vlgr_data_full_calendar[,m], na.rm=FALSE )
  vlgr_data_full_calendar_interpolate_4[,m] = zoo::na.spline( vlgr_data_full_calendar[,m], na.rm=FALSE, 1:nrow(vlgr_data_full_calendar) )
  print(m)
}

#' Check to see what interpolations look like
plot( vlgr_data_df[,1] , vlgr_data_df[,7], typ="l" , col="red")
lines( vlgr_data_full_calendar[,1] , vlgr_data_full_calendar_interpolate_1[,7], typ="l" , col="black")
lines( vlgr_data_full_calendar[,1] , vlgr_data_full_calendar_interpolate_2[,7], typ="l" , col="blue")
lines( vlgr_data_full_calendar[,1] , vlgr_data_full_calendar_interpolate_3[,7], typ="l" , col="green")
lines( vlgr_data_full_calendar[,1] , vlgr_data_full_calendar_interpolate_4[,7], typ="l" , col="violet")
#' Feels like it is better to use linear interpolation as the splines may be adding too much in some places

#' Repeat interpolation process with d(Rt)/dt data
#' Select dates with no d(Rt)/dt data and set data points as NA
missing_dates_drt_dt = subset(calendar_dates[,2],  !(calendar_dates[,2] %in% drt_dt[,"data"]) )
empty_mat_drt_dt = matrix( NA , nrow = length( missing_dates_drt_dt ) , ncol = 3 )
missing_dates_drt_dt = cbind( missing_dates_drt_dt , as.data.frame( empty_mat_drt_dt ) )
names(missing_dates_drt_dt) <- c("time","date","drt_dt","drt_dt_interpolate")
#' Merge the dates with and without d(Rt)/dt data

#drt_dt_temp = as.data.frame( drt_dt[,3] )
#drt_dt_temp = lubridate::decimal_date( as.numeric( drt_dt_temp ) )
#drt_dt_temp = cbind( drt_dt , lubridate::decimal_date(as.numeric(as.data.frame(drt_dt[,3]))) )
drt_dt_temp = as.data.frame( drt_dt )
drt_dt_temp[,1] =  drt_dt_temp[,3]
drt_dt_temp[,3] =  drt_dt_temp[,4]
drt_dt_temp[,c(2,4)] = NA
drt_dt_temp[,5:ncol(drt_dt_temp)]=NULL
names(drt_dt_temp) = names(missing_dates_drt_dt)

drt_dt_full_calendar = rbind( drt_dt_temp , missing_dates_drt_dt )
drt_dt_full_calendar <- data.table::setorder( drt_dt_full_calendar , cols = "time" )

#' Interpolate Var(LGR) for dates with no values
drt_dt_full_calendar[,4] = zoo::na.approx( drt_dt_full_calendar[,3], na.rm=FALSE )
drt_dt_full_calendar[,5] = zoo::na.approx( drt_dt_full_calendar[,3], na.rm=FALSE , 1:nrow(drt_dt_full_calendar) )
drt_dt_full_calendar[,6] = zoo::na.spline( drt_dt_full_calendar[,3], na.rm=FALSE )
drt_dt_full_calendar[,7] = zoo::na.spline( drt_dt_full_calendar[,3], na.rm=FALSE, 1:nrow(drt_dt_full_calendar) )
drt_dt_full_calendar[,8] = rowMeans( drt_dt_full_calendar[,c(4,5,6,7)] , na.rm = TRUE)

#' Check to see what interpolations look like
plot( drt_dt_full_calendar[,1] , drt_dt_full_calendar[,3], typ="p" , col="red")
lines( drt_dt_full_calendar[,1] , drt_dt_full_calendar[,4], typ="l" , col="black")
lines( drt_dt_full_calendar[,1] , drt_dt_full_calendar[,5], typ="l" , col="blue")
lines( drt_dt_full_calendar[,1] , drt_dt_full_calendar[,6], typ="l" , col="green")
lines( drt_dt_full_calendar[,1] , drt_dt_full_calendar[,7], typ="l" , col="violet")
lines( drt_dt_full_calendar[,1] , drt_dt_full_calendar[,8], typ="l" , col="blue")
abline(h=0)
#' Feels like it is better to use linear interpolation as the splines may be adding too much in some places

#' trim to date range of Var(LGR)
drt_dt_full_calendar_trimmed = subset( drt_dt_full_calendar ,  drt_dt_full_calendar[,1] %in% lubridate::decimal_date( vlgr_data_full_calendar[,1] ))
#' check that dates match
sum(drt_dt_full_calendar_trimmed[,1] == lubridate::decimal_date( vlgr_data_full_calendar[,1] ) ) == nrow(drt_dt_full_calendar_trimmed)

#'** calculate correlations **
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_12/analysis/correlation/2023_01/") 



#' 7 days rolling correlation
run_cor_07 = vlgr_data_full_calendar ; run_cor_07[,2:ncol(run_cor_07)]=NA
n = 7
for (i in n : nrow( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
  for (j in 2 : ncol( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
    cor_start = i - n +1
    cor_end = i
    #pre_cor = data.frame( vlgr_data_full_calendar[ cor_start:cor_end , j ] , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = data.frame( vlgr_data_full_calendar_interpolate_1[ cor_start:cor_end , j ] 
                          , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = subset( pre_cor , !is.na(pre_cor[,1]) & !is.na(pre_cor[,2]) & !is.infinite(pre_cor[,1]) & !is.infinite(pre_cor[,2]) ) #' Exclude dates with NA/NaN/Inf for either Var(lgr) or d(Rt)/dt (these are NA for reasons other than not having a tree for that date)
    run_cor_07[ i , j ] =  cor( pre_cor[,1] , pre_cor[,2] )
    #message("Row ",i,", vlgr ",j)
  }
  message("Row ",i,", vlgr ",j)
  if ( ( i %% 50 == 0 ) | ( i == 593 ) ){
    saveRDS( run_cor_07 , paste0("run_cor_07_",i,".rds"))
  }
}

#' 14 days rolling correlation
run_cor_14 = vlgr_data_full_calendar ; run_cor_14[,2:ncol(run_cor_14)]=NA
n = 14
for (i in n : nrow( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
  for (j in 2 : ncol( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
    cor_start = i - n +1
    cor_end = i
    #pre_cor = data.frame( vlgr_data_full_calendar[ cor_start:cor_end , j ] , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = data.frame( vlgr_data_full_calendar_interpolate_1[ cor_start:cor_end , j ] 
                          , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = subset( pre_cor , !is.na(pre_cor[,1]) & !is.na(pre_cor[,2]) & !is.infinite(pre_cor[,1]) & !is.infinite(pre_cor[,2]) ) #' Exclude dates with NA/NaN/Inf for either Var(lgr) or d(Rt)/dt (these are NA for reasons other than not having a tree for that date)
    run_cor_14[ i , j ] =  cor( pre_cor[,1] , pre_cor[,2] )
    #message("Row ",i,", vlgr ",j)
  }
  message("Row ",i,", vlgr ",j)
  if ( ( i %% 50 == 0 ) | ( i == 593 ) ){
    saveRDS( run_cor_14 , paste0("run_cor_14_",i,".rds"))
  }
}

#' 21 days rolling correlation
run_cor_21 = vlgr_data_full_calendar ; run_cor_21[,2:ncol(run_cor_21)]=NA
n = 21
for (i in n : nrow( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
  for (j in 2 : ncol( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
    cor_start = i - n +1
    cor_end = i
    #pre_cor = data.frame( vlgr_data_full_calendar[ cor_start:cor_end , j ] , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = data.frame( vlgr_data_full_calendar_interpolate_1[ cor_start:cor_end , j ] 
                          , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = subset( pre_cor , !is.na(pre_cor[,1]) & !is.na(pre_cor[,2]) & !is.infinite(pre_cor[,1]) & !is.infinite(pre_cor[,2]) ) #' Exclude dates with NA/NaN/Inf for either Var(lgr) or d(Rt)/dt (these are NA for reasons other than not having a tree for that date)
    run_cor_21[ i , j ] =  cor( pre_cor[,1] , pre_cor[,2] )
    #message("Row ",i,", vlgr ",j)
  }
  message("Row ",i,", vlgr ",j)
  if ( ( i %% 50 == 0 ) | ( i == 593 ) ){
    saveRDS( run_cor_21 , paste0("run_cor_21_",i,".rds"))
  }
}

#' 28 days rolling correlation
run_cor_28 = vlgr_data_full_calendar ; run_cor_28[,2:ncol(run_cor_28)]=NA
n = 28
for (i in n : nrow( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
  for (j in 2 : ncol( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
    cor_start = i - n +1
    cor_end = i
    #pre_cor = data.frame( vlgr_data_full_calendar[ cor_start:cor_end , j ] , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = data.frame( vlgr_data_full_calendar_interpolate_1[ cor_start:cor_end , j ] 
                          , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = subset( pre_cor , !is.na(pre_cor[,1]) & !is.na(pre_cor[,2]) & !is.infinite(pre_cor[,1]) & !is.infinite(pre_cor[,2]) ) #' Exclude dates with NA/NaN/Inf for either Var(lgr) or d(Rt)/dt (these are NA for reasons other than not having a tree for that date)
    run_cor_28[ i , j ] =  cor( pre_cor[,1] , pre_cor[,2] )
    #message("Row ",i,", vlgr ",j)
  }
  message("Row ",i,", vlgr ",j)
  if ( ( i %% 50 == 0 ) | ( i == 593 ) ){
    saveRDS( run_cor_28 , paste0("run_cor_28_",i,".rds"))
  }
}

#' 35 days rolling correlation
run_cor_35 = vlgr_data_full_calendar ; run_cor_35[,2:ncol(run_cor_35)]=NA
n = 35
for (i in n : nrow( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
  for (j in 2 : ncol( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
    cor_start = i - n +1
    cor_end = i
    #pre_cor = data.frame( vlgr_data_full_calendar[ cor_start:cor_end , j ] , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = data.frame( vlgr_data_full_calendar_interpolate_1[ cor_start:cor_end , j ] 
                          , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = subset( pre_cor , !is.na(pre_cor[,1]) & !is.na(pre_cor[,2]) & !is.infinite(pre_cor[,1]) & !is.infinite(pre_cor[,2]) ) #' Exclude dates with NA/NaN/Inf for either Var(lgr) or d(Rt)/dt (these are NA for reasons other than not having a tree for that date)
    run_cor_35[ i , j ] =  cor( pre_cor[,1] , pre_cor[,2] )
    #message("Row ",i,", vlgr ",j)
  }
  message("Row ",i,", vlgr ",j)
  if ( ( i %% 50 == 0 ) | ( i == 593 ) ){
    saveRDS( run_cor_35 , paste0("run_cor_35_",i,".rds"))
  }
}

#' 56 days rolling correlation
run_cor_56 = vlgr_data_full_calendar ; run_cor_56[,2:ncol(run_cor_56)]=NA
n = 56
for (i in n : nrow( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
  for (j in 2 : ncol( vlgr_data_full_calendar_interpolate_1 ) ){ #vlgr_data_full_calendar ) ){
    cor_start = i - n +1
    cor_end = i
    #pre_cor = data.frame( vlgr_data_full_calendar[ cor_start:cor_end , j ] , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = data.frame( vlgr_data_full_calendar_interpolate_1[ cor_start:cor_end , j ] 
                          , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
    pre_cor = subset( pre_cor , !is.na(pre_cor[,1]) & !is.na(pre_cor[,2]) & !is.infinite(pre_cor[,1]) & !is.infinite(pre_cor[,2]) ) #' Exclude dates with NA/NaN/Inf for either Var(lgr) or d(Rt)/dt (these are NA for reasons other than not having a tree for that date)
    run_cor_56[ i , j ] =  cor( pre_cor[,1] , pre_cor[,2] )
    #message("Row ",i,", vlgr ",j)
  }
  message("Row ",i,", vlgr ",j)
  if ( ( i %% 50 == 0 ) | ( i == 593 ) ){
    saveRDS( run_cor_56 , paste0("run_cor_56_",i,".rds"))
  }
}

#' alternative using lapply **NEEDS UPDATING TO REFLECT CHANGES TO ABOVE**
run_cor_28 = vlgr_data_df ; run_cor_28[,2:ncol(run_cor_28)]=NA
vlgr_dfc_list = as.list( vlgr_data_full_calendar[,2:ncol(vlgr_data_full_calendar)] )

run_cor_28 = lapply( seq_along( vlgr_dfc_list ) , function(i){
  n = 28
  #for (i in n:nrow( vlgr_data_full_calendar ) ){
  for (m in n:length( vlgr_dfc_list[[1]] ) ){
      cor_start = m - n +1
      cor_end = m
      pre_cor = data.frame( vlgr_dfc_list[[i]][ cor_start:cor_end ] , drt_dt_full_calendar_trimmed[ cor_start:cor_end , 4 ] )
      pre_cor = subset( pre_cor , !is.na(pre_cor[,1]) & !is.na(pre_cor[,2]) )
      run_cor_28[ m , i ] =  cor( pre_cor[,1] , pre_cor[,2] )
      message("Row ",m,", vlgr ",i)
    }
  #}
} )

#' plot example results
example_vlgr_col = which(colnames( run_cor_28 ) == "tfps_vlgr_simple_samp_false_001_mina07_maxa56_mdperc" )
run_cor_07_extract = cbind( run_cor_07[,1] , run_cor_07[, example_vlgr_col ] )
run_cor_14_extract = cbind( run_cor_14[,1] , run_cor_14[, example_vlgr_col ] )
run_cor_21_extract = cbind( run_cor_21[,1] , run_cor_21[, example_vlgr_col ] )
run_cor_28_extract = cbind( run_cor_28[,1] , run_cor_28[, example_vlgr_col ] )
run_cor_35_extract = cbind( run_cor_35[,1] , run_cor_35[, example_vlgr_col ] )
run_cor_56_extract = cbind( run_cor_56[,1] , run_cor_56[, example_vlgr_col ] )
#' Var(LGR)
plot(as.Date( vlgr_data_full_calendar_interpolate_1[,1],origin="1970-01-01"),vlgr_data_full_calendar_interpolate_1[,example_vlgr_col],typ="l",ylab="",lty=1,xlab="Date")
abline(h=0)
#' d(Rt)/dt
plot(as.Date( drt_dt_full_calendar_trimmed[,1],origin="1970-01-01"),drt_dt_full_calendar_trimmed[,4],typ="l",ylab="d(Rt)/dt",lty=1,xlab="Date")
abline(h=0)

#' Correlation
#' 7-day rolling correlation
plot(as.Date( run_cor_07_extract[,1],origin="1970-01-01"),run_cor_07_extract[,2],typ="l",lty=1,ylim=c(-1,1),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(as.Date( run_cor_07_extract[,1],origin="1970-01-01"),replicate( nrow(run_cor_07_extract) , mean(run_cor_07_extract[,2],na.rm=TRUE)) , lty=2,col="blue")
#' 14-day rolling correlation
plot(as.Date( run_cor_14_extract[,1],origin="1970-01-01"),run_cor_14_extract[,2],typ="l",lty=1,ylim=c(-1,1),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(as.Date( run_cor_14_extract[,1],origin="1970-01-01"),replicate( nrow(run_cor_14_extract) ,mean(run_cor_14_extract[,2],na.rm=TRUE)) , lty=2,col="blue")
#' 21-day rolling correlation
plot(as.Date( run_cor_21_extract[,1],origin="1970-01-01"),run_cor_21_extract[,2],typ="l",lty=1,ylim=c(-1,1),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(as.Date( run_cor_21_extract[,1],origin="1970-01-01"),replicate( nrow(run_cor_21_extract) , mean(run_cor_21_extract[,2],na.rm=TRUE)) , lty=2,col="blue")
#' 28-day rolling correlation
plot(as.Date( run_cor_28_extract[,1],origin="1970-01-01"),run_cor_28_extract[,2],typ="l",lty=1,ylim=c(-1,1),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(as.Date( run_cor_28_extract[,1],origin="1970-01-01"),replicate( nrow(run_cor_28_extract) , mean(run_cor_28_extract[,2],na.rm=TRUE)) , lty=2,col="blue")
#' 35-day rolling correlation
plot(as.Date( run_cor_35_extract[,1],origin="1970-01-01"),run_cor_35_extract[,2],typ="l",lty=1,ylim=c(-1,1),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(as.Date( run_cor_35_extract[,1],origin="1970-01-01"),replicate( nrow(run_cor_35_extract) , mean(run_cor_35_extract[,2],na.rm=TRUE)) , lty=2,col="blue")
#' 56-day rolling correlation
plot(as.Date( run_cor_56_extract[,1],origin="1970-01-01"),run_cor_56_extract[,2],typ="l",lty=1,ylim=c(-1,1),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(as.Date( run_cor_56_extract[,1],origin="1970-01-01"),replicate( nrow(run_cor_56_extract) , mean(run_cor_56_extract[,2],na.rm=TRUE)) , lty=2,col="blue")


#' plot all results
plot(run_cor_28[,1],run_cor_28[,2],typ="l")
abline(h=0)
for ( i in 3:ncol(run_cor_28) ){
  lines(run_cor_28[,1],run_cor_28[,i])
}

#' plot min, max, mean, median
#' 7-day rolling correlation
run_cor_07_min = apply(run_cor_07[,2:4321], 1, FUN = min, na.rm = TRUE)
run_cor_07_max = apply(run_cor_07[,2:4321], 1, FUN = max, na.rm = TRUE)
run_cor_07_mean = apply(run_cor_07[,2:4321], 1, FUN = mean, na.rm = TRUE,)
run_cor_07_median = apply(run_cor_07[,2:4321], 1, FUN = median, na.rm = TRUE)
run_cor_07_mean_mean = mean( run_cor_07_mean , na.rm = TRUE)

plot(run_cor_07[,1],run_cor_07_min,typ="l",lty=2,ylim=c(-1,1.5),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(run_cor_07[,1],run_cor_07_max,typ="l",lty=2)
lines(run_cor_07[,1],run_cor_07_mean,typ="l",lty=1,col="red")
lines(run_cor_07[,1],run_cor_07_median,typ="l",lty=1,col="blue")
lines(run_cor_07[,1],replicate( nrow(run_cor_07),run_cor_07_mean_mean),typ="l",lty=3,col="black")

legend("topright"
       ,legend=c("Daily range","Daily mean","Daily median",paste("Mean of daily means = ",round(run_cor_07_mean_mean,3)))
       ,lty=c(2,1,1,3)
       ,col=c("black","red","blue","black")
       , cex = 1.2)
#' 14-day rolling correlation
run_cor_14_min = apply(run_cor_14[,2:4321], 1, FUN = min, na.rm = TRUE)
run_cor_14_max = apply(run_cor_14[,2:4321], 1, FUN = max, na.rm = TRUE)
run_cor_14_mean = apply(run_cor_14[,2:4321], 1, FUN = mean, na.rm = TRUE,)
run_cor_14_median = apply(run_cor_14[,2:4321], 1, FUN = median, na.rm = TRUE)
run_cor_14_mean_mean = mean( run_cor_14_mean , na.rm = TRUE)

plot(run_cor_14[,1],run_cor_14_min,typ="l",lty=2,ylim=c(-1,1.5),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(run_cor_14[,1],run_cor_14_max,typ="l",lty=2)
lines(run_cor_14[,1],run_cor_14_mean,typ="l",lty=1,col="red")
lines(run_cor_14[,1],run_cor_14_median,typ="l",lty=1,col="blue")
lines(run_cor_14[,1],replicate( nrow(run_cor_14),run_cor_14_mean_mean),typ="l",lty=3,col="black")

legend("topright"
       ,legend=c("Daily range","Daily mean","Daily median",paste("Mean of daily means = ",round(run_cor_14_mean_mean,3)))
       ,lty=c(2,1,1,3)
       ,col=c("black","red","blue","black")
       , cex = 1.2)
#' 21-day rolling correlation
run_cor_21_min = apply(run_cor_21[,2:4321], 1, FUN = min, na.rm = TRUE)
run_cor_21_max = apply(run_cor_21[,2:4321], 1, FUN = max, na.rm = TRUE)
run_cor_21_mean = apply(run_cor_21[,2:4321], 1, FUN = mean, na.rm = TRUE,)
run_cor_21_median = apply(run_cor_21[,2:4321], 1, FUN = median, na.rm = TRUE)
run_cor_21_mean_mean = mean( run_cor_21_mean , na.rm = TRUE)

plot(run_cor_21[,1],run_cor_21_min,typ="l",lty=2,ylim=c(-1,1.5),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(run_cor_21[,1],run_cor_21_max,typ="l",lty=2)
lines(run_cor_21[,1],run_cor_21_mean,typ="l",lty=1,col="red")
lines(run_cor_21[,1],run_cor_21_median,typ="l",lty=1,col="blue")
lines(run_cor_21[,1],replicate( nrow(run_cor_21),run_cor_21_mean_mean),typ="l",lty=3,col="black")

legend("topright"
       ,legend=c("Daily range","Daily mean","Daily median",paste("Mean of daily means = ",round(run_cor_21_mean_mean,3)))
       ,lty=c(2,1,1,3)
       ,col=c("black","red","blue","black")
       , cex = 1.2)
#' 28-day rolling correlation
run_cor_28_min = apply(run_cor_28[,2:4321], 1, FUN = min, na.rm = TRUE)
run_cor_28_max = apply(run_cor_28[,2:4321], 1, FUN = max, na.rm = TRUE)
run_cor_28_mean = apply(run_cor_28[,2:4321], 1, FUN = mean, na.rm = TRUE,)
run_cor_28_median = apply(run_cor_28[,2:4321], 1, FUN = median, na.rm = TRUE)
run_cor_28_mean_mean = mean( run_cor_28_mean , na.rm = TRUE)

plot(run_cor_28[,1],run_cor_28_min,typ="l",lty=2,ylim=c(-1,1.5),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(run_cor_28[,1],run_cor_28_max,typ="l",lty=2)
lines(run_cor_28[,1],run_cor_28_mean,typ="l",lty=1,col="red")
lines(run_cor_28[,1],run_cor_28_median,typ="l",lty=1,col="blue")
lines(run_cor_28[,1],replicate( nrow(run_cor_28),run_cor_28_mean_mean),typ="l",lty=3,col="black")

legend("topright"
       ,legend=c("Daily range","Daily mean","Daily median",paste("Mean of daily means = ",round(run_cor_28_mean_mean,3)))
       ,lty=c(2,1,1,3)
       ,col=c("black","red","blue","black")
       , cex = 1.2)
#' 35-day rolling correlation
run_cor_35_min = apply(run_cor_35[,2:4321], 1, FUN = min, na.rm = TRUE)
run_cor_35_max = apply(run_cor_35[,2:4321], 1, FUN = max, na.rm = TRUE)
run_cor_35_mean = apply(run_cor_35[,2:4321], 1, FUN = mean, na.rm = TRUE,)
run_cor_35_median = apply(run_cor_35[,2:4321], 1, FUN = median, na.rm = TRUE)
run_cor_35_mean_mean = mean( run_cor_35_mean , na.rm = TRUE)

plot(run_cor_35[,1],run_cor_35_min,typ="l",lty=2,ylim=c(-1,1.5),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(run_cor_35[,1],run_cor_35_max,typ="l",lty=2)
lines(run_cor_35[,1],run_cor_35_mean,typ="l",lty=1,col="red")
lines(run_cor_35[,1],run_cor_35_median,typ="l",lty=1,col="blue")
lines(run_cor_35[,1],replicate( nrow(run_cor_35),run_cor_35_mean_mean),typ="l",lty=3,col="black")

legend("topright"
       ,legend=c("Daily range","Daily mean","Daily median",paste("Mean of daily means = ",round(run_cor_35_mean_mean,3)))
       ,lty=c(2,1,1,3)
       ,col=c("black","red","blue","black")
       , cex = 1.2)
#' 56-day rolling correlation
run_cor_56_min = apply(run_cor_56[,2:4321], 1, FUN = min, na.rm = TRUE)
run_cor_56_max = apply(run_cor_56[,2:4321], 1, FUN = max, na.rm = TRUE)
run_cor_56_mean = apply(run_cor_56[,2:4321], 1, FUN = mean, na.rm = TRUE,)
run_cor_56_median = apply(run_cor_56[,2:4321], 1, FUN = median, na.rm = TRUE)
run_cor_56_mean_mean = mean( run_cor_56_mean , na.rm = TRUE)

plot(run_cor_56[,1],run_cor_56_min,typ="l",lty=2,ylim=c(-1,1.5),ylab="Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
lines(run_cor_56[,1],run_cor_56_max,typ="l",lty=2)
lines(run_cor_56[,1],run_cor_56_mean,typ="l",lty=1,col="red")
lines(run_cor_56[,1],run_cor_56_median,typ="l",lty=1,col="blue")
lines(run_cor_56[,1],replicate( nrow(run_cor_56),run_cor_56_mean_mean),typ="l",lty=3,col="black")

legend("topright"
       ,legend=c("Daily range","Daily mean","Daily median",paste("Mean of daily means = ",round(run_cor_56_mean_mean,3)))
       ,lty=c(2,1,1,3)
       ,col=c("black","red","blue","black")
       , cex = 1.2)

#' Plot surface/heatmap of time shifted correlations
library(plotly)
#' 7-day rolling correlation
run_cor_07_matrix = as.matrix( run_cor_07[,2:4321] , nrow = nrow( run_cor_07 ) , ncol = ncol( run_cor_07[,2:4321] ) )
fig <- plot_ly(z = ~run_cor_07_matrix)
fig <- fig %>% add_surface()
fig <- fig %>% layout(title = 'Manually Specified Labels'
                      ,  xaxis = list(title = 'leading indicator parameter set')
                      ,  yaxis = list(title = 'Number of days Var(LGR) shifted relative to d(Rt)/dt')
                      #,  zaxis = list(title = 'Pearson correlation coefficient')
                      , legend = list(title=list(text='<b> Pearson correlation coefficient </b>')))
fig
plot_ly(z = run_cor_07_matrix , type = "heatmap", x = seq(1,593,1) , y = seq(1,4320,1))
#' 14-day rolling correlation
run_cor_14_matrix = as.matrix( run_cor_14[,2:4321] , nrow = nrow( run_cor_14 ) , ncol = ncol( run_cor_14[,2:4321] ) )
plot_ly(z = run_cor_14_matrix , type = "heatmap", x = seq(1,593,1) , y = seq(1,4320,1))
#' 21-day rolling correlation
run_cor_21_matrix = as.matrix( run_cor_21[,2:4321] , nrow = nrow( run_cor_21 ) , ncol = ncol( run_cor_21[,2:4321] ) )
plot_ly(z = run_cor_21_matrix , type = "heatmap", x = seq(1,593,1) , y = seq(1,4320,1))
#' 28-day rolling correlation
run_cor_28_matrix = as.matrix( run_cor_28[,2:4321] , nrow = nrow( run_cor_28 ) , ncol = ncol( run_cor_28[,2:4321] ) )
plot_ly(z = run_cor_28_matrix , type = "heatmap", x = seq(1,593,1) , y = seq(1,4320,1))
#' 35-day rolling correlation
run_cor_35_matrix = as.matrix( run_cor_35[,2:4321] , nrow = nrow( run_cor_35 ) , ncol = ncol( run_cor_35[,2:4321] ) )
plot_ly(z = run_cor_35_matrix , type = "heatmap", x = seq(1,593,1) , y = seq(1,4320,1))
#' 14-day rolling correlation
run_cor_56_matrix = as.matrix( run_cor_56[,2:4321] , nrow = nrow( run_cor_56 ) , ncol = ncol( run_cor_56[,2:4321] ) )
plot_ly(z = run_cor_56_matrix , type = "heatmap", x = seq(1,593,1) , y = seq(1,4320,1))

#' Find Var(LGR) and rolling period with highest mean correlation 
run_cor_07_vlgr_means = apply(run_cor_07[,2:4321], 2, FUN = mean, na.rm = TRUE)
plot(run_cor_07_vlgr_means,xlab="Different Var(LGR) leading indicators",ylab="Mean Pearson correlation coefficient",main="7-day rolling",ylim=c(-1,1),cex.axis=1.5,cex.lab=1.5,cex.main=1.5) ; abline(h=0)
run_cor_14_vlgr_means = apply(run_cor_14[,2:4321], 2, FUN = mean, na.rm = TRUE)
plot(run_cor_14_vlgr_means,xlab="Different Var(LGR) leading indicators",ylab="Mean Pearson correlation coefficient",main="14-day rolling",ylim=c(-1,1),cex.axis=1.5,cex.lab=1.5,cex.main=1.5) ; abline(h=0)
run_cor_21_vlgr_means = apply(run_cor_21[,2:4321], 2, FUN = mean, na.rm = TRUE)
plot(run_cor_21_vlgr_means,xlab="Different Var(LGR) leading indicators",ylab="Mean Pearson correlation coefficient",main="21-day rolling",ylim=c(-1,1),cex.axis=1.5,cex.lab=1.5,cex.main=1.5) ; abline(h=0)
run_cor_28_vlgr_means = apply(run_cor_28[,2:4321], 2, FUN = mean, na.rm = TRUE)
plot(run_cor_28_vlgr_means,xlab="Different Var(LGR) leading indicators",ylab="Mean Pearson correlation coefficient",main="28-day rolling",ylim=c(-1,1),cex.axis=1.5,cex.lab=1.5,cex.main=1.5) ; abline(h=0)
run_cor_35_vlgr_means = apply(run_cor_35[,2:4321], 2, FUN = mean, na.rm = TRUE)
plot(run_cor_35_vlgr_means,xlab="Different Var(LGR) leading indicators",ylab="Mean Pearson correlation coefficient",main="35-day rolling",ylim=c(-1,1),cex.axis=1.5,cex.lab=1.5,cex.main=1.5) ; abline(h=0)
run_cor_56_vlgr_means = apply(run_cor_56[,2:4321], 2, FUN = mean, na.rm = TRUE)
plot(run_cor_56_vlgr_means,xlab="Different Var(LGR) leading indicators",ylab="Mean Pearson correlation coefficient",main="56-day rolling",ylim=c(-1,1),cex.axis=1.5,cex.lab=1.5,cex.main=1.5) ; abline(h=0)

run_cor_07_vlgr_max = max( run_cor_07_vlgr_means , na.rm = TRUE ) ; run_cor_07_max_col =  which( run_cor_07_vlgr_means == run_cor_07_vlgr_max ) ; message(round(run_cor_07_vlgr_max,3)," ",names(run_cor_07_max_col)," ",run_cor_07_max_col)
run_cor_14_vlgr_max = max( run_cor_14_vlgr_means , na.rm = TRUE ) ; run_cor_14_max_col =  which( run_cor_14_vlgr_means == run_cor_14_vlgr_max ) ; message(round(run_cor_14_vlgr_max,3)," ",names(run_cor_14_max_col)," ",run_cor_14_max_col)
run_cor_21_vlgr_max = max( run_cor_21_vlgr_means , na.rm = TRUE ) ; run_cor_21_max_col =  which( run_cor_21_vlgr_means == run_cor_21_vlgr_max ) ; message(round(run_cor_21_vlgr_max,3)," ",names(run_cor_21_max_col)," ",run_cor_21_max_col)
run_cor_28_vlgr_max = max( run_cor_28_vlgr_means , na.rm = TRUE ) ; run_cor_28_max_col =  which( run_cor_28_vlgr_means == run_cor_28_vlgr_max ) ; message(round(run_cor_28_vlgr_max,3)," ",names(run_cor_28_max_col)," ",run_cor_28_max_col)
run_cor_35_vlgr_max = max( run_cor_35_vlgr_means , na.rm = TRUE ) ; run_cor_35_max_col =  which( run_cor_35_vlgr_means == run_cor_35_vlgr_max ) ; message(round(run_cor_35_vlgr_max,3)," ",names(run_cor_35_max_col)," ",run_cor_35_max_col)
run_cor_56_vlgr_max = max( run_cor_56_vlgr_means , na.rm = TRUE ) ; run_cor_56_max_col =  which( run_cor_56_vlgr_means == run_cor_56_vlgr_max ) ; message(round(run_cor_56_vlgr_max,3)," ",names(run_cor_56_max_col)," ",run_cor_56_max_col)

plot(run_cor_35[,1],run_cor_35[,1582+1],typ="l",lty=1,ylim=c(-1,1.5),ylab="Rolling Pearson correlation coefficient",xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)

#' Var(LGR)
plot(as.Date( vlgr_data_full_calendar_interpolate_1[,1],origin="1970-01-01"),vlgr_data_full_calendar_interpolate_1[,1582+1],typ="l",ylab="",lty=1,xlab="Date",cex.axis=1.5,cex.lab=1.5)
abline(h=0)

####    Save results
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_12/analysis/correlation") 
