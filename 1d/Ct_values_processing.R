#' Author: Kieran Drake
#' Date: June 2022
#' Purpose: Produce time series of cycle threshold values for PCR tests for three
#' SARS-CoV-2 gene targets

################################
# install required packages
if(!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") 
if(!requireNamespace("Rtools", quietly = TRUE))
  install.packages("Rtools")
if(!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
library(data.table)
if(!requireNamespace("R.utils", quietly = TRUE))
  install.packages('R.utils')
if(!requireNamespace("dplyr", quietly = TRUE))
  install.packages('dplyr')
if(!requireNamespace("asbio", quietly = TRUE))
  install.packages('asbio')
################################

# Set working directory where file is located
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values")

# See Methodology.doc for file source
sgtf_df <- fread("rtm_incoming_sgtf_20220609-185931-e31b9dff_sgtf_deaths.csv")
# sgtf_100_df <- sgtf_df[1:100,]

# Checking data
unique(sgtf_df$nhser_name) # Regions indicate that data is for England only
unique(sgtf_df$sgtf_under30ct) # results are: NA, 0, -9, 1

#' Cycle threshold data - see link below for details of terminology
#' https://rdrr.io/github/terminological/uk-covid-datatools/src/R/SPIMDatasetProvider.R
# p2ch1q: number of cycles for ORF1ab
# p2ch2q: number of cycles for N gene
# p2ch3q: number of cycles for S gene
# p2ch4q: number of cycles for control

# Plot histograms for each gene target and control
hist(sgtf_df$p2ch1cq,xlab = "ORF1ab Ct value")
hist(sgtf_df$p2ch2cq,xlab = "N Ct value")
hist(sgtf_df$p2ch3cq,xlab = "S Ct value")
hist(sgtf_df$p2ch4cq,xlab = "Control Ct value")

# Filter data frame
# Remove entries where PCR control is NA
sgtf_filtered_df = subset(sgtf_df, !is.na(sgtf_df$p2ch4cq))
# Select required columns
sgtf_filtered_df = data.frame("Date" = sgtf_filtered_df$specimen_date
                             , "O_Ct" = sgtf_filtered_df$p2ch1cq
                             , "N_Ct" = sgtf_filtered_df$p2ch2cq
                             , "S_Ct" = sgtf_filtered_df$p2ch3cq
                             , "Control_Ct" = sgtf_filtered_df$p2ch4cq
                             , "pillar" = sgtf_filtered_df$pillar)

# Checks
nrow(sgtf_df) # 18,771,690
nrow(sgtf_filtered_df) # 7,871,759
nrow(subset(sgtf_filtered_df,is.na(sgtf_filtered_df$O_Ct))) # 50,236
nrow(subset(sgtf_filtered_df,!is.na(sgtf_filtered_df$O_Ct))) # 7,821,523
nrow(subset(sgtf_filtered_df,is.na(sgtf_filtered_df$O_Ct))) + nrow(subset(sgtf_filtered_df,!is.na(sgtf_filtered_df$O_Ct))) # 7,871,759
nrow(subset(sgtf_filtered_df,is.na(sgtf_filtered_df$N_Ct))) # 27,797
nrow(subset(sgtf_filtered_df,!is.na(sgtf_filtered_df$N_Ct))) # 7,843,962
nrow(subset(sgtf_filtered_df,is.na(sgtf_filtered_df$N_Ct))) + nrow(subset(sgtf_filtered_df,!is.na(sgtf_filtered_df$N_Ct))) # 7,871,759
nrow(subset(sgtf_filtered_df,is.na(sgtf_filtered_df$S_Ct))) # 3,334,578
nrow(subset(sgtf_filtered_df,!is.na(sgtf_filtered_df$S_Ct))) # 4,537,181
nrow(subset(sgtf_filtered_df,is.na(sgtf_filtered_df$S_Ct))) + nrow(subset(sgtf_filtered_df,!is.na(sgtf_filtered_df$S_Ct))) # 7,871,759
nrow(subset(sgtf_filtered_df,is.na(sgtf_filtered_df$Control_Ct))) # 0
nrow(subset(sgtf_df,is.na(sgtf_df$p2ch4cq))) # 10,899,931
nrow(subset(sgtf_filtered_df,!is.na(sgtf_filtered_df$Control_Ct))) # 7,871,759
nrow(subset(sgtf_df,is.na(sgtf_df$p2ch4cq))) + nrow(sgtf_filtered_df) # 18,771,690
nrow(subset(sgtf_filtered_df,sgtf_filtered_df$Control_Ct == 0)) # 0
nrow(subset(sgtf_filtered_df,sgtf_filtered_df$Control_Ct != 0)) # 7,871,759

# Remove original dataframe to reduce memory usage
rm(sgtf_df)

# Split by pillar
unique(sgtf_filtered_df$pillar) # Pillar 1, Pillar 2 and NA
sgtf_pillar1_df = subset(sgtf_filtered_df, sgtf_filtered_df$pillar == "Pillar 1")
sgtf_pillar2_df = subset(sgtf_filtered_df, sgtf_filtered_df$pillar == "Pillar 2")
# Remove pillar column
sgtf_pillar1_df$pillar <- NULL
sgtf_pillar2_df$pillar <- NULL

# Check
nrow(sgtf_filtered_df) # 7,871,759
nrow(sgtf_pillar1_df) + nrow(sgtf_pillar2_df) + nrow(subset(sgtf_filtered_df,is.na(pillar))) # 7,871,759
rm(sgtf_filtered_df) # to reduce memory required

################################################################################
# Plot histograms of Ct value frequency by gene and pillar
par(mfrow = c(2,4))
hist(sgtf_pillar1_df$O_Ct, col='blue',xlim = c(0,50),ylim = c(0,0.2),freq = FALSE,xlab = "ORF1ab gene Ct value", main = paste("Pillar 1"))
hist(sgtf_pillar1_df$N_Ct, col='blue',xlim = c(0,50),ylim = c(0,0.2),freq = FALSE,xlab = "N gene Ct value",main=paste())
hist(sgtf_pillar1_df$S_Ct, col='blue',xlim = c(0,50),ylim = c(0,0.2),freq = FALSE,xlab = "S gene Ct value",main=paste())
hist(sgtf_pillar1_df$Control_Ct, col='blue',xlim = c(0,50),ylim = c(0,0.2),freq = FALSE,xlab = "Control Ct value",main=paste())

hist(sgtf_pillar2_df$O_Ct, col='green',xlim = c(0,50),ylim = c(0,0.2),freq = FALSE,xlab = "ORF1ab gene Ct value",main=paste("Pillar 2"))
hist(sgtf_pillar2_df$N_Ct, col='green',xlim = c(0,50),ylim = c(0,0.2),freq = FALSE,xlab = "N gene Ct value",main=paste())
hist(sgtf_pillar2_df$S_Ct, col='green',xlim = c(0,50),ylim = c(0,0.2),freq = FALSE,xlab = "S gene Ct value",main=paste())
hist(sgtf_pillar2_df$Control_Ct, col='green',xlim = c(0,50),ylim = c(0,0.2),freq = FALSE,xlab = "Control Ct value",main=paste())

# If want to add a second hist plot to the same axes use 'add=T'
par(mfrow=c(2,2))
# ORF1ab
hist(sgtf_pillar1_df$O_Ct
     , col=rgb(0,0,1,0.5) # red,green,blue and transparency
     , xlim = c(0,50)
     , ylim = c(0,0.2)
     , freq = FALSE
     , xlab = "Ct value"
     , main = paste("ORF1ab gene"))
hist(sgtf_pillar2_df$O_Ct, col=rgb(0,1,0,0.5), add=T, freq = FALSE)
legend("topright",legend=c("Pillar 1","Pillar 2"),col=c(rgb(0,0,1,0.5),rgb(0,1,0,0.5)),lty=1,lwd=5)
# N gene
hist(sgtf_pillar1_df$N_Ct
     , col=rgb(0,0,1,0.5) # red,green,blue and transparency
     , xlim = c(0,50)
     , ylim = c(0,0.2)
     , freq = FALSE
     , xlab = "Ct value"
     , main = paste("N gene"))
hist(sgtf_pillar2_df$N_Ct, col=rgb(0,1,0,0.5), add=T, freq = FALSE)
legend("topright",legend=c("Pillar 1","Pillar 2"),col=c(rgb(0,0,1,0.5),rgb(0,1,0,0.5)),lty=1,lwd=5)
# S gene
hist(sgtf_pillar1_df$S_Ct
     , col=rgb(0,0,1,0.5) # red,green,blue and transparency
     , xlim = c(0,50)
     , ylim = c(0,0.2)
     , freq = FALSE
     , xlab = "Ct value"
     , main = paste("S gene"))
hist(sgtf_pillar2_df$S_Ct, col=rgb(0,1,0,0.5), add=T, freq = FALSE)
legend("topright",legend=c("Pillar 1","Pillar 2"),col=c(rgb(0,0,1,0.5),rgb(0,1,0,0.5)),lty=1,lwd=5)
# Control
hist(sgtf_pillar1_df$Control_Ct
     , col=rgb(0,0,1,0.5) # red,green,blue and transparency
     , xlim = c(0,50)
     , ylim = c(0,0.2)
     , freq = FALSE
     , xlab = "Ct value"
     , main = paste("Control"))
hist(sgtf_pillar2_df$Control_Ct, col=rgb(0,1,0,0.5), add=T, freq = FALSE)
legend("topright",legend=c("Pillar 1","Pillar 2"),col=c(rgb(0,0,1,0.5),rgb(0,1,0,0.5)),lty=1,lwd=5)

################################################################################
# Plot tests by date
par(mfrow=c(1,1))
hist(sgtf_pillar1_df$Date
     , col=rgb(0,0,1,0.5) # red, green, blue and transparency
     #, xlim = c(0,50)
     , ylim = c(0,0.01)
     , breaks = 24
     , freq = FALSE
     , xlab = "Date"
     , main = paste("PCR tests with gene targets recorded"))
hist(sgtf_pillar2_df$Date, breaks = 24,col=rgb(0,1,0,0.5), add=T, freq = FALSE)
legend("topright",legend=c("Pillar 1","Pillar 2"),col=c(rgb(0,0,1,0.5),rgb(0,1,0,0.5)),lty=1,lwd=5)

################################################################################
#' Take the minimum (~highest viral load) and mean of the three gene target Ct values for
#' each test as some have gene dropout. Then average (mean and median) per day.
temp_df = sgtf_pillar1_df[2:4]
sgtf_pillar1_df$Ct_min <- apply(temp_df, 1, FUN = min, na.rm = TRUE)
sgtf_pillar1_df$Ct_mean <- apply(temp_df, 1, FUN = mean, na.rm = TRUE)
rm(temp_df)
temp_df = sgtf_pillar2_df[2:4]
sgtf_pillar2_df$Ct_min <- apply(temp_df, 1, FUN = min, na.rm = TRUE)
sgtf_pillar2_df$Ct_mean <- apply(temp_df, 1, FUN = mean, na.rm = TRUE)
rm(temp_df)
#' Where all three gene targets are NA then Ct_min = Inf and Ct_mean = NaN
#' so need to change these to NA before aggregating by date below
#' Pillar 1
sgtf_pillar1_df[sgtf_pillar1_df == Inf] <- NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
sgtf_pillar1_df[is.nan(sgtf_pillar1_df)] <- NA
# Pillar 2
sgtf_pillar2_df[sgtf_pillar2_df == Inf] <- NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
sgtf_pillar2_df[is.nan(sgtf_pillar2_df)] <- NA

# Calculate normalised (for Ct Control value) Ct values for three genes
# Pillar 1
sgtf_pillar1_df$O_Ct_norm <- sgtf_pillar1_df$O_Ct - sgtf_pillar1_df$Control_Ct
sgtf_pillar1_df$N_Ct_norm <- sgtf_pillar1_df$N_Ct - sgtf_pillar1_df$Control_Ct
sgtf_pillar1_df$S_Ct_norm <- sgtf_pillar1_df$S_Ct - sgtf_pillar1_df$Control_Ct
# Pillar 2
sgtf_pillar2_df$O_Ct_norm <- sgtf_pillar2_df$O_Ct - sgtf_pillar2_df$Control_Ct
sgtf_pillar2_df$N_Ct_norm <- sgtf_pillar2_df$N_Ct - sgtf_pillar2_df$Control_Ct
sgtf_pillar2_df$S_Ct_norm <- sgtf_pillar2_df$S_Ct - sgtf_pillar2_df$Control_Ct

# Calculate viral load proxy ln(-2^(Ct - Control)) 
# Dahdou et al (2020) https://doi.org/10.1016/j.jinf.2020.10.017
# Pillar 1
sgtf_pillar1_df$O_vl <- log( 2 ^ (-sgtf_pillar1_df$O_Ct_norm ))
sgtf_pillar1_df$N_vl <- log( 2 ^ (-sgtf_pillar1_df$N_Ct_norm ))
sgtf_pillar1_df$S_vl <- log( 2 ^ (-sgtf_pillar1_df$S_Ct_norm ))
# Pillar 2
sgtf_pillar2_df$O_vl <- log( 2 ^ (-sgtf_pillar2_df$O_Ct_norm ))
sgtf_pillar2_df$N_vl <- log( 2 ^ (-sgtf_pillar2_df$N_Ct_norm ))
sgtf_pillar2_df$S_vl <- log( 2 ^ (-sgtf_pillar2_df$S_Ct_norm ))

#' Take the maximum (~highest viral load) and mean of the three gene target viral load 
#' proxy values Ct values for each test as some have gene dropout. 
#' Then average (mean and median) per day.
temp_df = sgtf_pillar1_df[11:13]
sgtf_pillar1_df$vl_max <- apply(temp_df, 1, FUN = max, na.rm = TRUE)
sgtf_pillar1_df$vl_mean <- apply(temp_df, 1, FUN = mean, na.rm = TRUE)
rm(temp_df)
temp_df = sgtf_pillar2_df[11:13]
sgtf_pillar2_df$vl_max <- apply(temp_df, 1, FUN = max, na.rm = TRUE)
sgtf_pillar2_df$vl_mean <- apply(temp_df, 1, FUN = mean, na.rm = TRUE)
rm(temp_df)
#' Where all three gene targets are NA then vl_max = Inf and vl_mean = NaN
#' so need to change these to NA before aggregating by date below
#' Pillar 1
sgtf_pillar1_df[sgtf_pillar1_df == Inf] <- NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
sgtf_pillar1_df[is.nan(sgtf_pillar1_df)] <- NA
# Pillar 2
sgtf_pillar2_df[sgtf_pillar2_df == Inf] <- NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
sgtf_pillar2_df[is.nan(sgtf_pillar2_df)] <- NA


# Group by date, creating mean and median values
# Pillar 1 mean by date
Ct_p1_mean_df <- aggregate(sgtf_pillar1_df[,-1] 
                      , by = list(sgtf_pillar1_df$Date)
                      , FUN = "mean"
                      , na.rm = TRUE
                      )
colnames(Ct_p1_mean_df)[1] <- "Date"
# Pillar 2 mean by date
Ct_p2_mean_df <- aggregate(sgtf_pillar2_df[,-1] 
                      , by = list(sgtf_pillar2_df$Date)
                      , FUN = "mean"
                      , na.rm = TRUE
                      )
colnames(Ct_p2_mean_df)[1] <- "Date"
# Pillar 1 median by date
Ct_p1_median_df <- aggregate(sgtf_pillar1_df[,-1] #O_Ct + N_Ct + S_Ct + Control_Ct + Ct_min + Ct_mean ~ Date 
                           , by = list(sgtf_pillar1_df$Date)
                           , FUN = "median"
                           , na.rm = TRUE
                           )
colnames(Ct_p1_median_df)[1] <- "Date"
# Pillar 2 median by date
Ct_p2_median_df <- aggregate(sgtf_pillar2_df[,-1] #O_Ct + N_Ct + S_Ct + Control_Ct + Ct_min + Ct_mean ~ Date 
                           , by = list(sgtf_pillar2_df$Date)
                           , FUN = "median"
                           , na.rm = TRUE
)
colnames(Ct_p2_median_df)[1] <- "Date"

# plot mean
par(mfrow=c(1,1))
plot(Ct_p1_mean_df$Date,Ct_p1_mean_df$Ct_min,typ="l",col="blue")
lines(Ct_p1_mean_df$Date,Ct_p1_mean_df$Ct_mean,typ="l",col="darkgreen")
lines(Ct_p1_mean_df$Date,Ct_p1_mean_df$O_Ct,typ="l",col="lightblue")
lines(Ct_p1_mean_df$Date,Ct_p1_mean_df$N_Ct,typ="l",col="green")
lines(Ct_p1_mean_df$Date,Ct_p1_mean_df$S_Ct,typ="l",col="red")

plot(Ct_p2_mean_df$Date,Ct_p2_mean_df$Ct_min,typ="l",col="blue")
lines(Ct_p2_mean_df$Date,Ct_p2_mean_df$Ct_mean,typ="l",col="darkgreen")
lines(Ct_p2_mean_df$Date,Ct_p2_mean_df$O_Ct,typ="l",col="lightblue")
lines(Ct_p2_mean_df$Date,Ct_p2_mean_df$N_Ct,typ="l",col="green")
lines(Ct_p2_mean_df$Date,Ct_p2_mean_df$S_Ct,typ="l",col="red")

# plot median
par(mfrow=c(1,1))
plot(Ct_p1_median_df$Date,Ct_p1_median_df$Ct_min,typ="l",col="blue")
lines(Ct_p1_median_df$Date,Ct_p1_median_df$Ct_mean,typ="l",col="darkgreen")
lines(Ct_p1_median_df$Date,Ct_p1_median_df$O_Ct,typ="l",col="lightblue")
lines(Ct_p1_median_df$Date,Ct_p1_median_df$N_Ct,typ="l",col="green")
lines(Ct_p1_median_df$Date,Ct_p1_median_df$S_Ct,typ="l",col="red")

plot(Ct_p2_median_df$Date,Ct_p2_median_df$Ct_min,typ="l",col="blue")
lines(Ct_p2_median_df$Date,Ct_p2_median_df$Ct_mean,typ="l",col="darkgreen")
lines(Ct_p2_median_df$Date,Ct_p2_median_df$O_Ct,typ="l",col="lightblue")
lines(Ct_p2_median_df$Date,Ct_p2_median_df$N_Ct,typ="l",col="green")
lines(Ct_p2_median_df$Date,Ct_p2_median_df$S_Ct,typ="l",col="red")

# Difference between mean and min per day
x1 = Ct_p1_mean_df$Date
x2 = Ct_p2_mean_df$Date
y_p2_mean_diff = Ct_p2_mean_df$Ct_mean - Ct_p2_mean_df$Ct_min
y_p2_median_diff = Ct_p2_median_df$Ct_mean - Ct_p2_median_df$Ct_min

plot( x2 , y_p2_mean_diff , typ="l" , xlab = "Date")
lines( x2 , y_p2_median_diff , typ="l" , col="blue")

y_p2_median_mean_diff = y_p2_median_diff - y_p2_mean_diff
plot( x2 , y_p2_median_mean_diff , typ="l" , xlab = "Date" , ylim = c(-0.2, +0.3) )

plot(x1,Ct_p1_mean_df$Ct_mean - Ct_p1_mean_df$Ct_min,typ="l", col="red")
lines(x1,Ct_p1_median_df$Ct_mean - Ct_p1_median_df$Ct_min,typ="l",col="green")

# Load hospitalisations data, trim and organise
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs'
filename <- "data_2022-May-09 - UK Covid-19 cases by sample date.csv"
dat_type <- "cases"
#or
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
dat_type <- "hospitalisations"
cas_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)
#test
hosp_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)


plot( x2 , y_p2_median_mean_diff , typ="l" , xlab = "Date" , ylim = c(-0.2, +0.3) )
lines(hosp_df$date,hosp_df$cases/max(hosp_df$cases)*0.3,col="red")
lines(cas_df$date,cas_df$cases/max(cas_df$cases)*0.3,col="blue")


plot( x2 +60 , frollmean(y_p2_median_mean_diff, 21) , typ="l" , xlab = "Date" , ylim = c(-0.2, +0.3) )
lines(hosp_df$date,frollmean(hosp_df$cases/max(hosp_df$cases)*0.3, 1),col="red")
lines(cas_df$date,cas_df$cases/max(cas_df$cases)*0.3,col="blue")

par(mfrow=c(2,1))
plot(cas_df$Date,cases_df$cases/max(cases_df$cases),col="red")
plot(Ct_mean_mean_pillar2_df$Date
     ,Ct_mean_mean_pillar2_df$Ct_mean/max(Ct_mean_mean_pillar2_df$Ct_mean)
     ,typ="l"
     ,col="lightblue"
     ,ylim=c(0.55,0.65))

#' Calculate skew and standard deviation of daily distribution of patient min and
#' mean Ct value
#' Calc for Pillar 2, mean grouped by date
#' max viral load (and min Ct value)
Ct_skew_min_list = c()
Ct_st_dev_min_list = c()
vl_skew_max_list = c()
vl_st_dev_max_list = c()
for (i in 1:length( Ct_p2_mean_df$Date ) ) {
  print( i )
  Ct_s <- asbio::skew( subset( sgtf_pillar2_df , Date == Ct_p2_mean_df$Date[i] )$Ct_min )
  vl_s <- asbio::skew( subset( sgtf_pillar2_df , Date == Ct_p2_mean_df$Date[i] )$vl_max )
  Ct_st_dev <- sd( subset( sgtf_pillar2_df , Date == Ct_p2_mean_df$Date[i] )$Ct_min )
  vl_st_dev <- sd( subset( sgtf_pillar2_df , Date == Ct_p2_mean_df$Date[i] )$vl_max )
  Ct_skew_min_list = c( Ct_skew_min_list , Ct_s )
  vl_skew_max_list = c( vl_skew_max_list , vl_s )
  Ct_st_dev_min_list = c( Ct_st_dev_min_list , Ct_st_dev )
  vl_st_dev_max_list = c( vl_st_dev_max_list , vl_st_dev )
}
Ct_p2_mean_df$Ct_min_skew  = Ct_skew_min_list
Ct_p2_mean_df$Ct_min_stdev = Ct_st_dev_min_list
Ct_p2_mean_df$vl_max_skew  = vl_skew_max_list
Ct_p2_mean_df$vl_max_stdev = vl_st_dev_max_list

# mean
Ct_skew_mean_list = c()
Ct_st_dev_mean_list = c()
vl_skew_mean_list = c()
vl_st_dev_mean_list = c()
for (j in 1:length(Ct_p2_mean_df$Date)){
  print(j)
  Ct_s <- asbio::skew( subset( sgtf_pillar2_df , Date == Ct_p2_mean_df$Date[j] )$Ct_mean )
  Ct_st_dev <- sd( subset( sgtf_pillar2_df , Date == Ct_p2_mean_df$Date[j] )$Ct_mean )
  Ct_skew_mean_list = c( Ct_skew_mean_list , Ct_s )
  Ct_st_dev_mean_list = c( Ct_st_dev_mean_list , Ct_st_dev )
  vl_s <- asbio::skew( subset( sgtf_pillar2_df , Date == Ct_p2_mean_df$Date[j] )$vl_mean )
  vl_st_dev <- sd( subset( sgtf_pillar2_df , Date == Ct_p2_mean_df$Date[j] )$vl_mean )
  vl_skew_mean_list = c( vl_skew_mean_list , vl_s )
  vl_st_dev_mean_list = c( vl_st_dev_mean_list , vl_st_dev )
}
Ct_p2_mean_df$Ct_mean_skew = Ct_skew_mean_list
Ct_p2_mean_df$Ct_mean_stdev = Ct_st_dev_mean_list
Ct_p2_mean_df$vl_mean_skew = vl_skew_mean_list
Ct_p2_mean_df$vl_mean_stdev = vl_st_dev_mean_list

# Calc for Pillar 2, median grouped by date
#' max viral load (and min Ct value)
Ct_skew_min_list = c()
Ct_st_dev_min_list = c()
vl_skew_max_list = c()
vl_st_dev_max_list = c()
for (i in 1:length(Ct_p2_median_df$Date)){
  print(i)
  Ct_s <- asbio::skew( subset( sgtf_pillar2_df , Date == Ct_p2_median_df$Date[i] )$Ct_min )
  Ct_st_dev <- sd( subset( sgtf_pillar2_df , Date == Ct_p2_median_df$Date[i] )$Ct_min )
  Ct_skew_min_list = c( Ct_skew_min_list , Ct_s )
  Ct_st_dev_min_list = c( Ct_st_dev_min_list , Ct_st_dev )
  vl_s <- asbio::skew( subset( sgtf_pillar2_df , Date == Ct_p2_median_df$Date[i] )$vl_max )
  vl_st_dev <- sd( subset( sgtf_pillar2_df , Date == Ct_p2_median_df$Date[i] )$vl_max )
  vl_skew_max_list = c( vl_skew_max_list , vl_s )
  vl_st_dev_max_list = c( vl_st_dev_max_list , vl_st_dev )
}
Ct_p2_median_df$Ct_min_skew = Ct_skew_min_list
Ct_p2_median_df$Ct_min_stdev = Ct_st_dev_min_list
Ct_p2_median_df$vl_max_skew = vl_skew_max_list
Ct_p2_median_df$vl_max_stdev = vl_st_dev_max_list

Ct_skew_mean_list = c()
Ct_st_dev_mean_list = c()
vl_skew_mean_list = c()
vl_st_dev_mean_list = c()
for (j in 1:length(Ct_p2_mean_df$Date)){
  print(j)
  Ct_s <- asbio::skew( subset( sgtf_pillar2_df , Date == Ct_p2_median_df$Date[j] )$Ct_mean )
  Ct_st_dev <- sd( subset( sgtf_pillar2_df , Date == Ct_p2_median_df$Date[j] )$Ct_mean )
  Ct_skew_mean_list = c( Ct_skew_mean_list , Ct_s )
  Ct_st_dev_mean_list = c( Ct_st_dev_mean_list , Ct_st_dev )
  vl_s <- asbio::skew( subset( sgtf_pillar2_df , Date == Ct_p2_median_df$Date[j] )$vl_mean )
  vl_st_dev <- sd( subset( sgtf_pillar2_df,Date == Ct_p2_median_df$Date[j] )$vl_mean )
  vl_skew_mean_list = c( vl_skew_mean_list , vl_s )
  vl_st_dev_mean_list = c( vl_st_dev_mean_list , vl_st_dev )
}
Ct_p2_median_df$Ct_mean_skew = Ct_skew_mean_list
Ct_p2_median_df$Ct_mean_stdev = Ct_st_dev_mean_list
Ct_p2_median_df$vl_mean_skew = vl_skew_mean_list
Ct_p2_median_df$vl_mean_stdev = vl_st_dev_mean_list

rm(i , j
   , Ct_skew_mean_list , Ct_skew_min_list 
   , Ct_st_dev_mean_list , Ct_st_dev_min_list
   , vl_skew_mean_list , vl_skew_min_list 
   , vl_st_dev_mean_list , vl_st_dev_min_list
   , vl_skew_max_list , vl_st_dev_max_list )

# Write data to files
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files")
write.csv(sgtf_pillar1_df, file="C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files/sgtf_pillar1_df.csv")
write.csv(sgtf_pillar2_df, file="C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files/sgtf_pillar2_df_v2.csv")
write.csv(Ct_p1_mean_df, file="C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files/Ct_p1_mean_df.csv")
write.csv(Ct_p1_median_df, file="C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files/Ct_p1_median_df.csv")
write.csv(Ct_p2_mean_df, file="C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files/Ct_p2_mean_df_v2.csv")
write.csv(Ct_p2_median_df, file="C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/England Ct values/Processed files/Ct_p2_median_df_v2.csv")


# plot data
par(mfrow=c(1,1))
plot(Ct_p2_mean_df$Date,1-frollmean(Ct_p2_mean_df$min_skew,21),typ="l",ylim=c(0,1))
lines(Ct_p2_mean_df$Date,1-frollmean(Ct_p2_mean_df$min_stdev,21),col="green")
lines(hosp_df$date,frollmean(hosp_df$cases/max(hosp_df$cases), 21),col="red")
lines(cas_df$date,frollmean(cas_df$cases/max(cas_df$cases),21),col="blue")

plot(Ct_p2_mean_df$Date,1-frollmean(Ct_p2_mean_df$mean_stdev,21),typ="l",ylim=c(0,1))
lines(Ct_p2_mean_df$Date,1-frollmean(Ct_p2_mean_df$mean_skew,21),col="green")
lines(hosp_df$date,frollmean(hosp_df$cases/max(hosp_df$cases), 21),col="red")
lines(cas_df$date,frollmean(cas_df$cases/max(cas_df$cases),21),col="blue")

hist(subset(sgtf_pillar2_df,Date == "2020-12-20")$Ct_min,breaks= 12)
hist(subset(sgtf_pillar2_df,Date == "2021-03-20")$Ct_min,breaks= 12)
hist(subset(sgtf_pillar2_df,Date == "2021-12-20")$Ct_min,breaks= 12)
hist(subset(sgtf_pillar2_df,Date == "2022-04-24")$Ct_min,breaks= 12)

plot(Ct_p2_mean_df$min_skew,Ct_p2_mean_df$Ct_min)
plot(Ct_p2_mean_df$mean_skew,Ct_p2_mean_df$Ct_mean)
plot(Ct_p2_mean_df$O_Ct,Ct_p2_mean_df$N_Ct)
plot(Ct_p2_mean_df$O_Ct,Ct_p2_mean_df$S_Ct)
plot(Ct_p2_mean_df$N_Ct,Ct_p2_mean_df$S_Ct)