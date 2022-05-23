#' KD 18 May 2022
#' Simplified overlay
#' Load data for calculation of wave start/end dates
#' and manipulate results from growth methods

#' Uses non-standard functions, which need to be loaded: data_load(), 
#' gam_fitting(), wave_define_growth_method()
#' 
#' ###############################
install.packages("Rtools")
install.packages("data.table")
library(data.table)
install.packages('R.utils')
install.packages("TTR")
###############################

# Load data, trim and organise
folder <- 'C:/Users/kdrake/Documents/GitHub/Early_Warning_Signal'
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
dat_type <- "hospitalisations"
hosp_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)

#################################
# Create dataframe for storing summary information for different models
# !!!!WARNING!!!! This will remove any existing data in this dataframe
wave_define_df <- data.frame( "Data_type" = NA#c()
                              , "GAM method" = NA#c()
                              , "GAM smooth function" = NA#c()
                              , "Degrees_of_freedom_k" = NA#c()
                              , "k'" = NA#c()
                              , "k-index" = NA#c()
                              , "R_squared" = NA#c()
                              , "edf" = NA#c()
                              , "edf/k' ratio" = NA#c()
                              , "GCV" = NA#c()
                              , "Wave ID method" = NA#c()
                              , "Derivative_threshold" = NA#c()
                              , "Derivative_significance_level" = NA#c()
                              , "Growth_threshold" = NA#c()
                              , "Number_of_waves" = NA#c()
)
for (i in 1:100){
  j = i
  col_name = paste("wave",j,"start")
  additional_col = data.frame(placeholder_name = as.Date(NA))
  names(additional_col) <- col_name
  wave_define_df = cbind(wave_define_df,additional_col[1])
}
for (r in 1:100){
  s = r
  col_name = paste("wave",r,"end")
  additional_col = data.frame(placeholder_name = as.Date(NA))
  names(additional_col) <- col_name
  wave_define_df = cbind(wave_define_df,additional_col[1])
}

################################

# Define generalised additive model (GAM) inputs
GAM_smooth_function = "tp"
deg_free_k = 50 # Degrees of freedom for GAM model with regard to date

# Calculate generalised additive model (GAM)
cases_df = hosp_df
m <- gam_fitting(cases_df, GAM_smooth_function, deg_free_k)
# Display information for checking model quality and summary data
par(mfrow=c(2,2))
gam.check(m)
summary(m)
################################
# Return model from using growth of log(smoothed cases) method to identify wave dates
#growth_threshold = 0.030
counter = 0
for (gsf in c("tp","ts","ds","ps","cp","re","gp","cr","cs","cc","mrf")){
  for (k in seq(50,150,10)) {
    m <- gam_fitting(cases_df = hosp_df, GAM_smooth_function = gsf, deg_free_k = k)
    for (gt in seq(0,0.1,0.01)) {
      UK_model = growth_method(   m
                                , cases_df
                                , dat_type
                                , wave_bands_df
                                , GAM_smooth_function = gsf
                                , growth_threshold = gt
                                , deg_free_k = k
      )
      # Add output to dataframe collating previously produced model summaries
      wave_define_df = rbind(wave_define_df,UK_model[1,])
      counter = counter + 1
      message(counter," of 1331 complete")
    }
  }
}



# Remove the first row which contains NAs (from when first created)
# (Only needs to be done once)
wave_define_df <- wave_define_df[-c(1), ]

# Write dataframe to file
write.csv(wave_define_df_7waves, file="wave_define_hosp_growth_7waves_2022-05-23.csv")
##############################

# Plot data
par(mfrow=c(2,1))
par(mfrow=c(1,1))

# Plot cases and model and annotate with waves
plot(cases_df$time,cases_df$cases,xlab = "year",ylab = "cases",typ="l")
lines(cases_df$time,m$fitted.values,col="blue")
points(cases_df$time[wave_start_ix],m$fitted.values[wave_start_ix],col="red",cex=2)
#points(cases_df$time[wave_reset_ix],m$fitted.values[wave_reset_ix],col="green",cex=2)
legend("topright",c("cases","model","wave start","wave reset"),col=c("black","blue","red","darkgreen"),lty=1)
text(cases_df$time[wave_start_ix],m$fitted.values[wave_start_ix]+500,c(seq(1,length(wave_start_ix))),col="red")
#text(cases_df$time[wave_reset_ix],m$fitted.values[wave_reset_ix]+500,c(seq(1,length(wave_reset_ix))),col="darkgreen")

# Plot residuals
plot(cases_df$time,m$residuals,xlab = "year", ylab = "GAM residuals")

# Plot growth and log(growth)
plot(growth_df$date,growth_df$growth,xlab="date",ylab="growth for smoothed log cases")
lines(growth_df$date,y=replicate(length(growth_df$date),0),col="black")
points(growth_df$date[wave_start_ix],growth_df$growth[wave_start_ix],col="red")
plot(cases_df$time,cases_df$smoothedLogCases,xlab="date",ylab="smoothed log cases")
lines(cases_df$time,y=replicate(length(cases_df$time),0),col="black")
points(cases_df$time[wave_start_ix],cases_df$smoothedLogCases[wave_start_ix],col="red")

# Plot smoothed log(cases) and growth in smoothed log(cases) and mark on wave dates
par(mfrow=c(1,1))
plot(cases_df$date,cases_df$smoothedLogCases,xlab="date",ylab="") #,ylab="smoothed log cases")
lines(cases_df$date,y=replicate(length(cases_df$date),0),col="black")
line_for_xaxis = replicate(length(cases_df$date),0)
points(cases_df$date[wave_start_ix],line_for_xaxis[wave_start_ix],col="red")
points(cases_df$date[wave_reset_ix],line_for_xaxis[wave_reset_ix],col="darkgreen")
#points(cases_df$date[wave_start_ix],cases_df$smoothedLogCases[wave_start_ix],col="red")
lines(growth_df$date,growth_df$growth*20,col="blue") # growth in smoothed log cases
lines(growth_df$date,growth_df$growth_adj*20,col="blue") # growth in smoothed log cases
legend("bottomright",legend = c("Smoothed log(cases)","growth in smoothed log(cases) * 20","Wave start","Wave reset"),col=c("black","blue","red","darkgreen"),lty=1,cex=0.8) 
text(as.Date(cases_df$date[wave_start_ix]),growth_df$growth[wave_start_ix]*20+0.5,c(seq(1,length(wave_start_ix))),col="red")
text(as.Date(cases_df$date[wave_reset_ix]),growth_df$growth[wave_reset_ix]*20+0.5,c(seq(1,length(wave_reset_ix))),col="darkgreen")
lines(growth_df$date,replicate(length(growth_df$growth),growth_threshold*20),col="blue")

#growth_rate = TTR::ROC(m$fitted.values, type = "discrete")
#plot(date,growth_rate)

# Just looking at models that idenitified 7 wave starts
# Plots 
par(mfrow=c(3,2))
# wave 1
hist(wave_define_df_7waves$'wave 1 start', breaks=50,xlab = "Date")
plot(wave_define_df_7waves$'wave 1 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 1 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 1 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 1 start',wave_define_df_7waves$R_squared,col="green")

# wave 2
hist(wave_define_df_7waves$'wave 2 start', breaks=50,xlab = "Date")
plot(wave_define_df_7waves$'wave 2 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 2 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 2 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 2 start',wave_define_df_7waves$R_squared,col="green")

# wave 3
hist(wave_define_df_7waves$'wave 3 start', breaks=50,xlab = "Date")
plot(wave_define_df_7waves$'wave 3 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 3 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 3 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 3 start',wave_define_df_7waves$R_squared,col="green")

# wave 4
hist(wave_define_df_7waves$'wave 4 start', breaks=50,xlab = "Date")
plot(wave_define_df_7waves$'wave 4 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 4 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 4 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 4 start',wave_define_df_7waves$R_squared,col="green")

# wave 5
hist(wave_define_df_7waves$'wave 5 start', breaks=50,xlab = "Date")
plot(wave_define_df_7waves$'wave 5 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 5 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 5 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 5 start',wave_define_df_7waves$R_squared,col="green")

# wave 6
hist(wave_define_df_7waves$'wave 6 start', breaks=50,xlab = "Date")
plot(wave_define_df_7waves$'wave 6 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 6 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 6 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 6 start',wave_define_df_7waves$R_squared,col="green")

# wave 7
hist(wave_define_df_7waves$'wave 7 start', breaks=50,xlab = "Date")
plot(wave_define_df_7waves$'wave 7 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 7 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 7 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 7 start',wave_define_df_7waves$R_squared,col="green")

par(mfrow=c(5,7))
# histogram of dates
hist(wave_define_df_7waves$'wave 1 start', breaks=50)
hist(wave_define_df_7waves$'wave 2 start', breaks=50)
hist(wave_define_df_7waves$'wave 3 start', breaks=50)
hist(wave_define_df_7waves$'wave 4 start', breaks=50)
hist(wave_define_df_7waves$'wave 5 start', breaks=50)
hist(wave_define_df_7waves$'wave 6 start', breaks=50)
hist(wave_define_df_7waves$'wave 7 start', breaks=50)
# edf/k ratio
plot(wave_define_df_7waves$'wave 1 start',wave_define_df_7waves$edf.k..ratio,col="blue",ylab = "edf/k ratio")
plot(wave_define_df_7waves$'wave 2 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 3 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 4 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 5 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 6 start',wave_define_df_7waves$edf.k..ratio,col="blue")
plot(wave_define_df_7waves$'wave 7 start',wave_define_df_7waves$edf.k..ratio,col="blue")
#growth threshold
plot(wave_define_df_7waves$'wave 1 start',wave_define_df_7waves$Growth_threshold,col="purple",ylab = "Growth threshold")
plot(wave_define_df_7waves$'wave 2 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 3 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 4 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 5 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 6 start',wave_define_df_7waves$Growth_threshold,col="purple")
plot(wave_define_df_7waves$'wave 7 start',wave_define_df_7waves$Growth_threshold,col="purple")
# k-index
plot(wave_define_df_7waves$'wave 1 start',wave_define_df_7waves$k.index,col="red",ylab="k-index")
plot(wave_define_df_7waves$'wave 2 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 3 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 4 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 5 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 6 start',wave_define_df_7waves$k.index,col="red")
plot(wave_define_df_7waves$'wave 7 start',wave_define_df_7waves$k.index,col="red")
# R^2
plot(wave_define_df_7waves$'wave 1 start',wave_define_df_7waves$R_squared,col="green",ylab="R^2")
plot(wave_define_df_7waves$'wave 2 start',wave_define_df_7waves$R_squared,col="green")
plot(wave_define_df_7waves$'wave 3 start',wave_define_df_7waves$R_squared,col="green")
plot(wave_define_df_7waves$'wave 4 start',wave_define_df_7waves$R_squared,col="green")
plot(wave_define_df_7waves$'wave 5 start',wave_define_df_7waves$R_squared,col="green")
plot(wave_define_df_7waves$'wave 6 start',wave_define_df_7waves$R_squared,col="green")
plot(wave_define_df_7waves$'wave 7 start',wave_define_df_7waves$R_squared,col="green")