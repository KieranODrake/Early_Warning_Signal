# Author: Kieran Drake
# Date: May 2022
# Produce time series of positive and negative pillar 2 PCR tests for Covid-19
# and calculate positivity rate time series from line list data

#########################
# install packages for fread (which is quicker than read_csv)
if(!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") 
if(!requireNamespace("Rtools", quietly = TRUE))
  install.packages("Rtools")
if(!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
if(!requireNamespace("R.utils", quietly = TRUE))
  install.packages('R.utils')
if(!requireNamespace("dplyr", quietly = TRUE))
  install.packages('dplyr')

library(BiocManager)
library(data.table)
library(R.utils)
library(dplyr)

#########################
# Set working directory where files located
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/Positivity rates/Raw/")

# Case data from line list
positive_df <- fread("Anonymised Combined Line List 20220523.csv") #https://imperiallondon.sharepoint.com/sites/ncov/Shared%20Documents/Forms/AllItems.aspx?e=5%3Ac1d1fb4fb3bd473dad6051fee991f62b&at=9&RootFolder=%2Fsites%2Fncov%2FShared%20Documents%2F2019%2DnCoV%2FData%20collection%2FUK%2FLinelist&FolderCTID=0x012000D3E93DEBFB0FF44E8C182382B05245D2
negative_df <- fread("20220417 TOTAL Negative tests.csv") # https://imperiallondon.sharepoint.com/sites/ncov/Shared%20Documents/Forms/AllItems.aspx?e=5%3A24082fdfe2274ccf8e0e843e34cf94c7&at=9&RootFolder=%2Fsites%2Fncov%2FShared%20Documents%2F2019%2DnCoV%2FData%20collection%2FUK%2FNegatives%5FPillars&FolderCTID=0x012000D3E93DEBFB0FF44E8C182382B05245D2
View( positive_df[ 1 : 100 , ] )
View( negative_df[ 1 : 100 , ] )

# Trim for Pillar 2
positive_df <- subset( positive_df , positive_df$pillar == "PILLAR2" )
negative_df <- subset( negative_df , negative_df$pillar == 2 )

# Trim non-PCR
positive_df <- subset( positive_df , positive_df$case_category %in% c( "PCR_ONLY" , "LFT_WITHPCT" ) )
negative_df <- subset( negative_df , is.na( negative_df$LFT_Flag ) )

# Select columns
positive_df <- data.frame(   "specimen_date" = as.Date( positive_df$specimen_date , format = "%d/%m/%Y" ) 
                           , "NHSER_name"    = positive_df$NHSER_name 
                           )
negative_df <- data.frame(  "earliestspecimendate" = negative_df$earliestspecimendate
                          , "NHSRegion"            = negative_df$NHSRegion
                          , "total_negative"       = as.numeric( negative_df$total_negative )
                          )

# Check regions
unique( positive_df$NHSER_name ) # Only regions in England returned plus 17024 with blank NHSER_name
unique( negative_df$NHSRegion ) # Only regions in England returned

# reshape from long to wide
#negative_df$total_negative <- as.numeric(negative_df$total_negative)
negativetest_df <- aggregate( . ~ NHSRegion + earliestspecimendate , negative_df , sum )
negativetest2_df <- reshape( negativetest_df , idvar = "earliestspecimendate" , timevar = "NHSRegion" , direction = "wide" )
#mutate(negative_p2_wide_df$Eng_total = rowSums (.))

# test count of negatives retained after reshaping
sum( negative_df$total_negative ) # 100383019
sum( negativetest_df$total_negative , na.rm = TRUE ) # 100383019
sum( negativetest2_df[ , 2 : 8 ] , na.rm = TRUE ) # 100383019
negative_df = negativetest2_df
rm( negativetest_df )
rm( negativetest2_df )

# Reshape and test count of positives retained after reshaping # 11,437,743 rows before
positivetest_df = data.frame( rbind( table( positive_df ) ) )
sum( positivetest_df , na.rm = TRUE ) # 11,437,743 after
positive_df = positivetest_df
rm( positivetest_df )
positive_df$Date = rownames( positive_df )

# Sum totals for England
negative_df$Eng_total <- rowSums( Filter( is.numeric , negative_df ) , na.rm = TRUE )
negative_out_df <- data.frame(   "Date" = as.Date( negative_df$earliestspecimendate ) 
                               , "Eng_total_negative" = negative_df$Eng_total 
                               )
positive_df$Eng_total <- rowSums( Filter( is.numeric , positive_df[ , 1:8 ] ) , na.rm = TRUE )
positive_out_df <- data.frame(  "Date" = as.Date( positive_df$Date )
                              , "Eng_total_positive" = positive_df$Eng_total 
                              )

# Merge positive and negative
library(dplyr)
positivity_rate_df = dplyr::left_join( positive_out_df , negative_out_df, by = "Date" )

# Calculate positivity rate
positivity_rate_df$Positivity_rate_England = positivity_rate_df$Eng_total_positive / ( positivity_rate_df$Eng_total_positive + positivity_rate_df$Eng_total_negative )
#' This version returns the mean for 7 values (3 before, 1 on, and 3 after the date) #' positivity_rate_df$pos_rate_7d_roll_avg_Eng = zoo::rollmean( positivity_rate_df$Positivity_rate_England , k = 7 , fill = NA )
positivity_rate_df$pos_rate_7d_roll_avg_Eng = data.table::frollmean( positivity_rate_df$Positivity_rate_England , n = 7 , fill = NA )

# Plot positivity rate
plot(  positivity_rate_df$Date
     , positivity_rate_df$Positivity_rate_England
     , typ="l"
     , ylim = c(0,0.5)
     , ylab = "Covid-19 PCR pillar 2 positivity rate - England"
     , xlab = " Date"
     , col = "red")
lines( positivity_rate_df$Date , positivity_rate_df$pos_rate_7d_roll_avg_Eng , typ = "l" , col = "black" )

# Write data to file
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs/Positivity rates/Outputs/")
write.csv( positivity_rate_df , file = "sc2_positivity_rate_England.csv" )
saveRDS( positivity_rate_df , file = "sc2_positivity_rate_England.rds" )

# Plot against hospitalisations - data loaded from 'wave define - overlay.R' and 'data_load.R'
#' Load case/hospitalisation data
folder <- 'C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Inputs'
#filename <- "data_2022-May-09 - UK Covid-19 cases by sample date.csv"
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
#dat_type <- "cases"
dat_type <- "hospitalisations"
#' data_load() @ C:\Users\kdrake\GitHub\Early_Warning_Signal
dat_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)
#'plot
plot( dat_df$date , dat_df$cases , xlab = "Date" , ylab = dat_type , typ="l" )
mult_factor = 15000 #* max(cases_df$cases)
lines(positivity_rate_df$Date
      , positivity_rate_df$Positivity_rate_England * mult_factor
      , typ="l"
      , col="red")
# add moving average to plot
#moving_avg = zoo::rollmean(positivity_rate_df$Positivity_rate_England, k = 7, fill = NA)
moving_avg = data.table::frollmean( positivity_rate_df$Positivity_rate_England , n = 7 , fill = NA )
lines(positivity_rate_df$Date
      , moving_avg * mult_factor
      , typ="l"
      , col="blue")
legend(c("UK hospitalisations"))

rm( negative_df , negative_out_df , positive_df , positive_out_df)
