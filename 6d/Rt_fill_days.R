#' Finding the dates at which Rt goes above 1 (critical transition)
#' 
#' ####    Load Rt data    ####
#' Rt data downloaded from https://www.gov.uk/guidance/the-r-value-and-growth-rate 221207_R_and_growth_rate_time_series_for_publication_v1.0.ods
#' Asterisks removed from some dates in date column and column reformatted to 'Date' in Excel
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_12/analysis/correlation") 
rt = read.csv("Rt.csv")
rt$mid <- ( rt$lower + rt$upper ) / 2
rt$date <- as.Date( rt$date , format = "%d/%m/%Y") # Change format of date
rt$time <- lubridate::decimal_date( rt$date ) # Add column for decimal date

#' This is a bi-weekly time series but we would like to know the exact date so
#' we fill in the intervening days
deg_free_k = 100
GAM_smooth_function = "tp"
m_mid <- mgcv::gam( mid ~ s( time, bs=GAM_smooth_function, k = deg_free_k), data = rt )
m_lower <- mgcv::gam( lower ~ s( time, bs=GAM_smooth_function, k = deg_free_k), data = rt )
m_upper <- mgcv::gam( upper ~ s( time, bs=GAM_smooth_function, k = deg_free_k), data = rt )

#' Plot
par(mfrow=c(2,1))
plot(rt$time,rt$mid,typ="l",lwd=2,col="blue",xlab="year",ylab="Rt",ylim=c(0.5,2.0))
lines( rt$time,rt$lower , typ="l" , col = "blue" , lwd=1, lty=2 )
lines( rt$time,rt$upper , typ="l" , col = "blue" , lwd=1, lty=2 )
lines( m_mid$model[[2]] , m_mid$fitted.values , typ="l" , col = "red" , lwd=2 )
legend("topright",legend=c("Rt (midpoint of 90% CI for England)", "Rt 90% CI band",paste0("Rt model fit (spline=",GAM_smooth_function,", k=",deg_free_k,")"),"Rt = 1")
       , lty=c(1,2,1,2), lwd=c(2,1,2,2),col=c("blue","blue","red","black")
)
abline(h=1,lty=2)
#' plot GAMs
plot( m_mid$model[[2]] , m_mid$fitted.values , typ="l",lwd=2,col="blue",xlab="year",ylab="Rt",ylim=c(0.5,2.0))
lines( m_lower$model[[2]] , m_lower$fitted.values , typ="l" , col = "blue" , lwd=1, lty=2 )
lines( m_upper$model[[2]] , m_upper$fitted.values , typ="l" , col = "blue" , lwd=1, lty=2 )
legend("topright",legend=c(paste0("Rt model fit (spline=",GAM_smooth_function,", k=",deg_free_k,")"), "Rt 90% CI band","Rt = 1")
       , lty=c(1,2,2), lwd=c(2,1,2,2),col=c("blue","blue","black")
)
abline(h=1,lty=2)

#' Full list of dates between start and end of bi-weekly time series
calendar_dates <- seq( min( rt$date ) , max( rt$date ) , 1 )
rt_extended <- data.frame(   "date"      = calendar_dates
                           , "lower"     = replicate( length( calendar_dates) , NA )
                           , "upper"     = replicate( length( calendar_dates) , NA )
                           , "mid"       = replicate( length( calendar_dates) , NA )
                           , "time"      = lubridate::decimal_date( calendar_dates )
                           , "GAM_lower" = replicate( length( calendar_dates) , NA )
                           , "GAM_upper" = replicate( length( calendar_dates) , NA )
                           , "GAM_mid"   = replicate( length( calendar_dates) , NA )
                          )
#' Only need dates not already in the data
'%ni%' = Negate( '%in%' )
rt_extended = subset( rt_extended , rt_extended$date %ni% rt$date )
#' Add GAM modelled values
rt_model = cbind( rt , m_lower$fitted.values, m_upper$fitted.values , m_mid$fitted.values  )
names( rt_model ) <- names( rt_extended )
rt_extended = rbind( rt_extended , rt_model )
rt_extended <- data.table::setorder( rt_extended , cols = "time" )
#' Then use xts::na.approx to fill in the empty dates
#' GAM lower CI
rt_extended[,9]  <- na.approx( rt_extended[,6] ) ; names(rt_extended)[9] <- "GAM_lower_approx"
#rt_extended[,10] <- na.approx(dates_extended[,2] , 1:400 )
rt_extended[,10] <- na.spline( rt_extended[,6] ) ; names(rt_extended)[10] <- "GAM_lower_spline"
#rt_extended[,6] <- na.spline(dates_extended[,2] , 1:400 )
#rt_extended[,7] <- rowMeans(dates_extended[,c(3,4,5,6)])
#' GAM upper CI
rt_extended[,11] <- na.approx( rt_extended[,7] ) ; names(rt_extended)[11] <- "GAM_upper_approx"
rt_extended[,12] <- na.spline( rt_extended[,7] ) ; names(rt_extended)[12] <- "GAM_upper_spline"
#' GAM upper CI
rt_extended[,13] <- na.approx( rt_extended[,8] ) ; names(rt_extended)[13] <- "GAM_mid_approx"
rt_extended[,14] <- na.spline( rt_extended[,8] ) ; names(rt_extended)[14] <- "GAM_mid_spline"

setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2022_12/analysis/correlation")
saveRDS( rt_extended , "Rt_data_extended.rds")
write.csv( rt_extended , "Rt_data_extended.csv")


#' Plot
par(mfrow=c(2,1))
plot( rt$time , rt$mid , typ="l" , lwd=2 , col="blue" , xlab = "year" , ylab = "Rt" , ylim = c( 0.5 , 2.0 ) )
lines( rt_extended$time , rt_extended$GAM_mid_approx , typ="l" , col = "red" , lwd=1, lty=2 )
lines( rt_extended$time , rt_extended$GAM_mid_spline , typ="l" , col = "green" , lwd=1, lty=2 )

lines( rt$time,rt$upper , typ="l" , col = "blue" , lwd=1, lty=2 )
lines( m_mid$model[[2]] , m_mid$fitted.values , typ="l" , col = "red" , lwd=2 )
legend("topright",legend=c("Rt (midpoint of 90% CI for England)", "Rt 90% CI band",paste0("Rt model fit (spline=",GAM_smooth_function,", k=",deg_free_k,")"),"Rt = 1")
       , lty=c(1,2,1,2), lwd=c(2,1,2,2),col=c("blue","blue","red","black")
)
abline(h=1,lty=2)
#' plot GAMs
plot( m_mid$model[[2]] , m_mid$fitted.values , typ="l",lwd=2,col="blue",xlab="year",ylab="Rt",ylim=c(0.5,2.0))
lines( m_lower$model[[2]] , m_lower$fitted.values , typ="l" , col = "blue" , lwd=1, lty=2 )
lines( m_upper$model[[2]] , m_upper$fitted.values , typ="l" , col = "blue" , lwd=1, lty=2 )
legend("topright",legend=c(paste0("Rt model fit (spline=",GAM_smooth_function,", k=",deg_free_k,")"), "Rt 90% CI band","Rt = 1")
       , lty=c(1,2,2), lwd=c(2,1,2,2),col=c("blue","blue","black")
)
abline(h=1,lty=2)


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