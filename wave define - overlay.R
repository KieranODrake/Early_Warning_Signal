#' KD 18 May 2022
#' Load data for calculation of wave start/end dates
#' and manipulate results from growth and derivative methods

#' Uses non-standard functions, which need to be loaded: data_load(), 
#' gam_fitting(), wave_define_derivative_method(), wave_define_growth_method()

###############################
# install packages for fread (which is quicker than read_csv)
install.packages("Rtools")
install.packages("data.table")
library(data.table)
install.packages('R.utils')
install.packages("TTR")
install.packages("gratia") #https://gavinsimpson.github.io/gratia/
###############################

# Load data, trim and organise
folder <- 'C:/Users/kdrake/Documents/GitHub/Early_Warning_Signal'

filename <- "data_2022-May-09 - UK Covid-19 cases by sample date.csv"
dat_type <- "cases"
#or
filename <- "data_2022-May-09 - UK Covid-19 hospital admissions.csv"
dat_type <- "hospitalisations"

cas_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)
#test
hosp_df <- data_load( folder , filename , dat_type ) # Outputs cases_df(date,cases,time,wday)
plot(cas_df$date,cas_df$cases)
lines(hosp_df$date,hosp_df$cases*10)

################################

# Define generalised additive model (GAM) inputs
GAM_smooth_function = "cc"
deg_free_k = 55 # Degrees of freedom for GAM model with regard to date

# Calculate generalised additive model (GAM)
cases_df = cas_df
cases_df = hosp_df
m <- gam_fitting(cases_df, GAM_smooth_function, deg_free_k)
# Display information for checking model quality and summary data
par(mfrow=c(2,2))
gam.check(m)
summary(m)

################################

# Define inputs for derivative and growth methods for defining wave start/end dates

# Manually define wave date bands 
# (semi-arbitrary earliest possible start date for wave i.e. just after previous wave peak) 
wave_bands_start <- c(as.Date(min(cases_df$date)),    # Wuhan
                      as.Date("2020-04-23"), # Alpha
                      as.Date("2020-11-14"), # Alpha
                      as.Date("2021-01-06"), # Delta
                      as.Date("2021-07-18"), # Delta
                      as.Date("2021-09-05"), # Delta
                      as.Date("2021-10-22"), # Omicron
                      as.Date("2022-01-05"), # Omicron
                      as.Date("2022-01-27"), # Omicron
                      as.Date("2022-03-25")  # Omicron
)
# (semi-arbitrary earliest possible end date for wave i.e. mid-point of next wave peak)
wave_bands_end <- c(as.Date("2020-03-27"), # Wuhan
                    as.Date("2020-10-10"),    # Wuhan
                    as.Date("2020-12-18"), # Alpha
                    as.Date("2021-06-29"), # Alpha
                    as.Date("2021-08-21"), # Delta
                    as.Date("2021-10-11"), # Delta
                    as.Date("2021-12-22"), # Delta
                    as.Date("2022-01-18"), # Omicron
                    as.Date("2022-03-07"), # Omicron
                    as.Date(max(cases_df$date))  # Omicron
)
wave_bands_df <- data.frame("start"= wave_bands_start,
                            "end" = wave_bands_end)
rm(wave_bands_start,wave_bands_end)

# Plot to show wave bands
par(mfrow=c(1,1))
plot(hosp_df$date,hosp_df$cases,xlim=c(as.Date("2020-01-01"),max(cases_df$date)),xlab="Date",ylab="Hospitalisations and Cases/40")
lines(cases_df$date,cases_df$cases/40,col="blue")
abline(v=wave_bands_df$start,col="red")
abline(v=wave_bands_df$end,col="green")
text(wave_bands_df$start+10,4000,seq(1,length(wave_bands_df$start)),col="red")
text(wave_bands_df$end-10,4200,seq(1,length(wave_bands_df$end)),col="green")


# Return model from using derivative method O'Brien & Clements (2021) and supplementary methods
derivative_threshold = 0
significance_level = 0.95 # Used to set the upper and lower bands in the for the derivative used to determine the wave phase
UK_model_derivative = derivative_method( m, cases_df, dat_type, wave_bands_df, GAM_smooth_function, threshold, significance_level, deg_free_k)

# Return model from using growth of log(smoothed cases) method
growth_threshold = 0
UK_model_growth = growth_method( m, cases_df, dat_type, wave_bands_df, GAM_smooth_function, threshold, significance_level, deg_free_k)

#################################

# Run function
start_time = Sys.time()

# Remove old lists from previous runs (if any)
rm(deg_free_k_list)
rm(R_squared_list)
rm(GCV_list)
rm(number_of_waves_list)
rm(wave_start_dates_list)
rm(method_list)
rm(wave_define)
rm(wave_define_df)

# (Re)Initialise lists to record data from function runs
deg_free_k_list = c()
R_squared_list = c()
GCV_list <- c()
number_of_waves_list = c()
wave_start_dates_list <- vector("list",46)
method_list = c()
wave_define <- vector("list",46)

# F statistic
# p-values

###########################################

# Create dataset using multiple k values in GAM
# !!!!WARNING!!!! This will remove any existing data in this dataframe
wave_define_df <- data.frame( "Data_type" = c()
                            , "GAM method" = c()
                            , "GAM smooth function" = c()
                            , "Degrees_of_freedom_k" = c()
                            , "k'" = c()
                            , "k-index" = c()
                            , "R_squared" = c()
                            , "edf" = c()
                            , "edf/k' ratio" = c()
                            , "GCV" = c()
                            , "Wave ID method" = c()
                            , "Derivative_threshold" = c()
                            , "Derivative_significance_level" = c()
                            , "Growth_threshold" = c()
                            , "Number_of_waves" = c()
                            , "Wave 1 start date" = c()
                            , "Wave 2 start date" = c()
                            , "Wave 3 start date" = c()
                            , "Wave 4 start date" = c()
                            , "Wave 5 start date" = c()
                            , "Wave 6 start date" = c()
                            , "Wave 7 start date" = c()
                            , "Wave 8 start date" = c()
                            , "Wave 9 start date" = c()
                            , "Wave 1 end date" = c()
                            , "Wave 2 end date" = c()
                            , "Wave 3 end date" = c()
                            , "Wave 4 end date" = c()
                            , "Wave 5 end date" = c()
                            , "Wave 6 end date" = c()
                            , "Wave 7 end date" = c()
                            , "Wave 8 end date" = c()
                            , "Wave 9 end date" = c()
)

# Create dataset using multiple k values in GAM
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

# Loop through range of k values for GAM
# Record R2, F and p values and wave start dates for each model

#######################

# set method
#method = "growth" 
#method = "derivative"

# Run function
k_index = 999 # there is an issue with k_index within the script - can check if 999 comes through into output dataframe
start_time = Sys.time()
#for (dat_type in c("cases","hospitalisations")){
  #dat_type = "hosp"
  #cases_df = data.frame(ifelse( dat_type == "cases" , cas_df[,] , hosp_df[,]))

  dat_type = "cases"
  cases_df = cas_df
  for (gsf in c("cc","tp","ts","ds","cr","cs")){ # already done "re","gp","ps","cp", and "so","sos" didn't work
    for (k in 5:300) {
      
      deg_free_k = k
      GAM_smooth_function = gsf
      
      # Calculate generalised additive model (GAM)
      m <- gam_fitting(cases_df, GAM_smooth_function, deg_free_k)
      
      ## Only look at models with a k-index > 1
      #gam_check_capture <- capture.output(gam.check(m))
      #k_index <- as.numeric(strsplit(gam_check_capture[13], " ")[[1]][7],5)
      #if (k_index > 1){
      
      
      #if (method == "growth"){
        # Identify waves using growth of log(smoothed cases) method
        UK_model = growth_method(   m
                                  , cases_df
                                  , dat_type
                                  , wave_bands_df
                                  , GAM_smooth_function
                                  , growth_threshold = 0
                                  , deg_free_k
                                )
      #}  
      
      # Add row to results dataframe
      wave_define_df = rbind(wave_define_df,UK_model[1,])
      
      #if (method == "derivative"){
        # Identify waves using derivative method
        UK_model = derivative_method( m
                                      , cases_df
                                      , dat_type
                                      , wave_bands_df
                                      , GAM_smooth_function
                                      , derivative_threshold = 0
                                      , significance_level = 0.95
                                      , deg_free_k
                                    )
      #}
        # Add row to results dataframe
        wave_define_df = rbind(wave_define_df,UK_model[1,])
    
        # message to track progress in for loop
        message("Data type: ", dat_type, "Smooth method: ", gsf, "k= ",k," k-index= ",k_index," Number of waves: ",UK_model$Number_of_waves)
        current_time = Sys.time()
        time_elapsed = as.numeric(current_time) - as.numeric(start_time)
        total_expected_time = 26400
        perc_complete = (time_elapsed/total_expected_time) * 100
        time_remaining = (total_expected_time - time_elapsed)/60
        message( "Percentage complete: ",perc_complete, " %. "," Approximate time remaining: ",time_remaining," minutes.")
        
      #} # for if statement above relating to k-index
    }
  }
#}

  dat_type = "hospitalisations"
  cases_df = hosp_df
  for (gsf in c("re","gp","cc","tp","ts","ds","cr","cs")){ # already done "ps","cp", and "so","sos" didn't work
    for (k in 5:300) {
      
      deg_free_k = k
      GAM_smooth_function = gsf
      
      # Calculate generalised additive model (GAM)
      m <- gam_fitting(cases_df, GAM_smooth_function, deg_free_k)
      
      ## Only look at models with a k-index > 1
      #gam_check_capture <- capture.output(gam.check(m))
      #k_index <- as.numeric(strsplit(gam_check_capture[13], " ")[[1]][7],5)
      #if (k_index > 1){
      
      
      #if (method == "growth"){
      # Identify waves using growth of log(smoothed cases) method
      UK_model = growth_method(   m
                                  , cases_df
                                  , dat_type
                                  , wave_bands_df
                                  , GAM_smooth_function
                                  , growth_threshold = 0
                                  , deg_free_k
      )
      #}  
      
      # Add row to results dataframe
      wave_define_df = rbind(wave_define_df,UK_model[1,])
      
      #if (method == "derivative"){
      # Identify waves using derivative method
      UK_model = derivative_method( m
                                    , cases_df
                                    , dat_type
                                    , wave_bands_df
                                    , GAM_smooth_function
                                    , derivative_threshold = 0
                                    , significance_level = 0.95
                                    , deg_free_k
      )
      #}
      # Add row to results dataframe
      wave_define_df = rbind(wave_define_df,UK_model[1,])
      
      # message to track progress in for loop
      message("Data type: ", dat_type, "Smooth method: ", gsf, "k= ",k," k-index= ",k_index," Number of waves: ",UK_model$Number_of_waves)
      current_time = Sys.time()
      time_elapsed = as.numeric(current_time) - as.numeric(start_time)
      total_expected_time = 26400
      perc_complete = (time_elapsed/total_expected_time) * 100
      time_remaining = (total_expected_time - time_elapsed)
      message( "Percentage complete: ",round(perc_complete,1), " %. "," Approximate time remaining: ",round(time_remaining,1)," minutes.")
      
      #} # for if statement above relating to k-index
    }
  }
  #}

end_time = Sys.time()
message("Total time taken: ", end_time - start_time)

#Remove the first row which if NAs
wave_define_df <- wave_define_df[-c(1), ]

# Convert dates
for (i in 16:16){
  for (j in 1:2){
    test[[i]][j] = ifelse( is.na(test[[i]][j])
                                                , NA
                                                , as.Date(test[[i]][j], origin = "1970-01-01")
                                              )
  }
}

for (i in 16:215){
  #for (j in 1:5928){
    test$Degrees_of_freedom_k[[i]] = as.Date(test[[i]], origin = "1970-01-01")
  #}
}

library(dplyr)


write.csv(wave_define_df, file="wave_define_df_2022-05-23.csv")

############################
# Plotting wave start dates

par(mfrow=c(1,1))
plot(wave_define_df$Wave.1.start.date,
     wave_define_df$Degrees_of_freedom_k,
     xlim=c(as.Date("2020-01-01"),as.Date("2022-04-30")),
     #ylim = c(0,600),
     xlab = "Date",
     ylab = "Degrees of freedom (k) in GAM")
points(wave_define_df$Wave.2.start.date,wave_define_df$Degrees_of_freedom_k,cex=0.5,col="grey")
points(wave_define_df$Wave.3.start.date,wave_define_df$Degrees_of_freedom_k,cex=0.5)
points(wave_define_df$Wave.4.start.date,wave_define_df$Degrees_of_freedom_k,cex=0.5,col="grey")
points(wave_define_df$Wave.5.start.date,wave_define_df$Degrees_of_freedom_k,cex=0.5)
points(wave_define_df$Wave.6.start.date,wave_define_df$Degrees_of_freedom_k,cex=0.5,col="grey")
points(wave_define_df$Wave.7.start.date,wave_define_df$Degrees_of_freedom_k,cex=0.5)
points(wave_define_df$Wave.8.start.date,wave_define_df$Degrees_of_freedom_k,cex=0.5,col="grey")
points(wave_define_df$Wave.9.start.date,wave_define_df$Degrees_of_freedom_k,cex=0.5)
lines(cases_df$date,cases_df$cases*(max(wave_define_df$Degrees_of_freedom_k)/max(cases_df$cases)),col="blue")


mean_date_decimal <- mean(lubridate::decimal_date(wave_define_df$Wave.9.start.date[!is.na(wave_define_df$Wave.9.start.date)]))
mean_date <- format(date_decimal(mean_date_decimal),"%Y-%m-%d")
y = replicate(length(wave_define_df$Degrees_of_freedom_k),mean_date)
abline(wave_define_df$Degrees_of_freedom_k,y,col="red")

# K against number of waves
derivative_only_df = wave_define_df[wave_define_df$Wave.ID.method == "Derivative", ] #subset(wave_define_df, Wave.ID.method="Derivative")
growth_only_df = subset(wave_define_df, Wave.ID.method=="Growth in log(smoothed cases)")
derivative_only_thresh_df = wave_define_df[(wave_define_df$Wave.ID.method == "Derivative") & (wave_define_df$Derivative_threshold == 0.1), ] #subset(wave_define_df, Wave.ID.method="Derivative")
growth_only_thresh_df = wave_define_df[(wave_define_df$Wave.ID.method == "Growth in log(smoothed cases)") & (wave_define_df$Growth_threshold == 0.01), ]

# plot all derivative and all growth
plot(derivative_only_df$Degrees_of_freedom_k,derivative_only_df$Number_of_waves,ylab="Number of waves identified",xlab="Degrees of freedom (k) in GAM")
points(growth_only_df$Degrees_of_freedom_k,growth_only_df$Number_of_waves,col="green")
legend("bottomright",legend=c("growth of log(smoothed cases) method","derivative method"),col=c("black","green"),lty=1)
# plot filtered threshold level
par(mfrow=c(1,1))
plot(derivative_only_df$Degrees_of_freedom_k,derivative_only_df$Number_of_waves,ylab="Number of waves identified",xlab="Degrees of freedom (k) in GAM")
points(growth_only_thresh_df$Degrees_of_freedom_k,growth_only_thresh_df$Number_of_waves,col="red")
points(growth_only_df$Degrees_of_freedom_k,growth_only_df$Number_of_waves,col="green")
points(derivative_only_thresh_df$Degrees_of_freedom_k,derivative_only_thresh_df$Number_of_waves,col="blue")

legend("bottomright",legend=c("growth of log(smoothed cases) method","derivative method"),col=c("black","green"),lty=1)



# K against number of R^2, k-index, edf/k ratio
plot(wave_define_df$Degrees_of_freedom_k, wave_define_df$k.index,xlab="Degrees of freedom (k) in GAM",ylab="")
points(wave_define_df$Degrees_of_freedom_k, wave_define_df$R_squared,xlab="Degrees of freedom (k) in GAM",col="blue")
points(wave_define_df$Degrees_of_freedom_k, wave_define_df$edf.k..ratio,xlab="Degrees of freedom (k) in GAM",col="red")
legend("bottomright",legend=c("R^2","k-index","edf / k ratio"),col=c("blue","black","red"),lty=1)


############################

# Plot summary data of variation of k = degrees of freedom in GAM
plot(x = wave_define_df$Degrees_of_freedom_in_GAM,
     y = wave_define_df$R_squared,
     xlim=c(0,50),
     ylim=c(0,1.2),
     xlab="Degrees of freedom in GAM",
     ylab = "")
points(x = wave_define_df$Degrees_of_freedom_in_GAM,
       y = 1-wave_define_df$GCV/max(wave_define_df$GCV),
       col="red")
lines(x = wave_define_df$Degrees_of_freedom_in_GAM,
      y = wave_define_df$number_of_waves/max(wave_define_df$number_of_waves),
      col="blue")
legend("bottomright",
       legend = c("R^2","1-GCV(max value normalised to 1)","No. of waves"),
       col=c("black","red","blue"),
       lty=1,
       cex = 1)
text(x = wave_define_df$Degrees_of_freedom_in_GAM,
     y = wave_define_df$number_of_waves/max(wave_define_df$number_of_waves)+0.03,
     labels = wave_define_df$number_of_waves,
     cex=1,
     col="blue")

################################

# OLD - PROBABLY NOT USEFUL
for (i in 1:length(wave_define)){
  #wave_define_df$Degrees_of_freedom_in_GAM[i] = wave_define[[i]][1]
  #wave_define_df$R_squared[i] = wave_define[[i]][2]
  #wave_define_df$GCV[i] = wave_define[[i]][3]
  #wave_define_df$number_of_waves[i] = wave_define[[i]][4]
  #wave_define_df$wave_start_dates[i] = wave_define[[i]][5]
  #wave_define_df$wave_identification_method[i] = wave_define[[i]][6]
  wave_define_df[nrow(wave_define_df) + 1] = c(wave_define[[i]][1],
                                               wave_define[[i]][2],
                                               wave_define[[i]][3],
                                               wave_define[[i]][4],
                                               #t(wave_define[[i]][5]),
                                               wave_define[[i]][6])
}

View(wave_define_df)
################################

# Plotting outputs

################################

# Define chart and data labels depending on data type loaded into function
y_axis_title <- ifelse(dat_type == "cases",
                       "Covid-19 cases in the UK",
                       "Covid-19 hospitalisations in the UK"
)

################################

# Plot outputs from derivative method

#current
#par(mfrow=c(1,1))
#plot(waves$year[which(!is.na(waves$cases_growing))],which(!is.na(waves$cases_growing)),col="red",type = "l" )

par(mfrow=c(2,2))
# Plot outputs
# Data and model
plot(waves$year,cases_df$cases,col="red",typ="p",xlab = "Year",ylab=y_axis_title)
lines(waves$year,m$fitted.values,col="blue",typ="l")
legend(legend = c("Data","Model"),col=c("red","blue"),lty=c(1,1),"topleft")
text(3.5,1250,labels = merge("R^2 = ",summary(m)$r.sq),col="black") 

par(mfrow=c(1,1))
# Plot cases by wave
g <- zoo(waves$cases_growing,waves$year) # zoo used so can plot discontinuous data
s <- zoo(waves$cases_stationary,waves$year)
d <- zoo(waves$cases_declining,waves$year)
plot(g,col="red",type = "l",xlab="Year",ylab= y_axis_title )
points(waves$year[wave_start_ix],cases_df$cases[wave_start_ix],col="red")
lines(s,col="black",type = "l" )
lines(d,col="green",type = "l" )
legend(legend = c("Positive","No trend","Negative"),col=c("red","black","green"),lty=1,title = "Predicted non-linear case trend","topleft") #,xjust =1,yjust=1)
text(cases_df$time[wave_start_ix],cases_df$cases[wave_start_ix]+500,c(seq(1,length(wave_start_ix))),col="red")
text(cases_df$time[wave_start_ix],cases_df$cases[wave_start_ix]+1000,labels = format(as.Date(wave_dates_nondecimal$start), "%d-%b"),cex=0.6,col="red")

# Plot residuals
plot(cases_df$time,m$residuals,xlab = "year", ylab = "GAM residuals")

library(gridExtra)
wd <-tableGrob(wave_dates_nondecimal)
grid.arrange(plot1,wd)
plot(wd)

# Derivative calculated from GAM model
plot(waves$year,waves$derivative_of_GAM,col="Blue",typ="l",xlab="Year",ylab="First derivative of GAM model")
lines(waves$year,y=replicate(length(waves$year),0))
lines(waves$year,waves$derivative_of_GAM_upper,lty=2,col="blue")
lines(waves$year,waves$derivative_of_GAM_lower,lty=2,col="blue")
points(waves$year[wave_start_ix],waves$derivative_of_GAM_lower[wave_start_ix],col="red")
conf_legend = paste(toString(significance_level*100),"% confidence band")
legend("bottomright", legend = c("First derivative",conf_legend,"Wave start"),col=c("blue","blue","red"),lty=c(1,2))
text(cases_df$time[wave_start_ix],waves$derivative_of_GAM_lower[wave_start_ix]+10,c(seq(1,length(wave_start_ix))),col="red")
text(waves$year[wave_start_ix],waves$derivative_of_GAM_lower[wave_start_ix]+5,labels = format(as.Date(wave_dates_nondecimal$start), "%d-%b"),cex=0.6,col="red") #labels = format(as.Date(wave_dates_nondecimal$start), "%d-%b-%Y")
# Plot original GAM derivative to check that derivative lines up with cases slope
#par(mfrow=c(1,1))
#plot(cases_df$time,cases_df$cases,ylab = "cases",xlab = "date",type="l")
#legend("top", legend = c("cases","GAM derivative"),col=c("black","blue"),lty = 1)
#lines(cases_df$time,replicate(length(GAM_derivative$derivative),50000))
#lines(cases_df$time,GAM_derivative$derivative*2000+50000, col = "blue")

############################

# Plot outputs from growth method

library(gridExtra)
wd <-tableGrob(wave_dates_nondecimal)
#grid.arrange(plot1,wd)
plot(wd)

# Plot data
par(mfrow=c(2,2))
par(mfrow=c(2,1))
par(mfrow=c(1,1))

# Plot cases and model and annotate with waves
plot(cases_df$time,cases_df$cases,xlab = "year",ylab = "cases",typ="l")
lines(cases_df$time,m$fitted.values,col="blue")
points(wave_dates$start,m$fitted.values[wave_start_ix],col="red") #, lab = wave_dates$start)
points(wave_dates$reset,m$fitted.values[wave_reset_ix],col="green")
legend("topright",c("cases","model","wave start","wave reset"),col=c("black","blue","red","darkgreen"),lty=1)
text(wave_dates$start,m$fitted.values[wave_start_ix]+500,c(seq(1,length(wave_start_ix))),col="red")
text(wave_dates$reset,m$fitted.values[wave_reset_ix]+500,c(seq(1,length(wave_reset_ix))),col="darkgreen")

# Plot residuals
plot(cases_df$time,m$residuals,xlab = "year", ylab = "GAM residuals")

# Plot growth and log(growth)
plot(growth_df$date,growth_df$growth,xlab="date",ylab="growth for smoothed log cases")
lines(growth_df$date,y=replicate(length(growth_df$date),0),col="black")
points(growth_df$date[wave_start_ix],growth_df$growth[wave_start_ix],col="red")
plot(cases_df$time,cases_df$smoothedLogCases,xlab="date",ylab="smoothed log cases")
lines(cases_df$time,y=replicate(length(cases_df$time),0),col="black")
points(cases_df$time[wave_start_ix],cases_df$smoothedLogCases[wave_start_ix],col="red")

# Plot cases (or hospitalisations) with model and highlight start of waves
par(mfrow=c(1,1))
plot(cases_df$time,cases_df$cases,xlab = "year",ylab = "cases",typ="l")
lines(cases_df$time,m$fitted.values,col="blue")
points(wave_dates$start,m$fitted.values[wave_start_ix],col="red") #, lab = wave_dates$start)
points(wave_dates$reset,m$fitted.values[wave_reset_ix],col="green") #, lab = wave_dates$start)
legend("top",legend = c("Data","Model","Wave start", "Wave reset"),col=c("black","blue","red","green"),lty=1) 
text(100,100,labels = merge("R^2 = ",summary(m)$r.sq),col="black") 

# Plot smoothed log(cases) and growth in smoothed log(cases) and mark on wave dates
par(mfrow=c(1,1))
plot(cases_df$date,cases_df$smoothedLogCases,xlab="date",ylab="") #,ylab="smoothed log cases")
lines(cases_df$date,y=replicate(length(cases_df$date),0),col="black")
line_for_xaxis = replicate(length(cases_df$date),0)
points(cases_df$date[wave_start_ix],line_for_xaxis[wave_start_ix],col="red")
points(cases_df$date[wave_reset_ix],line_for_xaxis[wave_reset_ix],col="darkgreen")
#points(cases_df$date[wave_start_ix],cases_df$smoothedLogCases[wave_start_ix],col="red")
lines(growth_df$date,growth_df$growth*20,col="blue") # growth in smoothed log cases
legend("bottomright",legend = c("Smoothed log(cases)","growth in smoothed log(cases) * 20","Wave start","Wave reset"),col=c("black","blue","red","darkgreen"),lty=1,cex=0.8) 
text(as.Date(wave_dates_nondecimal$start),growth_df$growth[wave_start_ix]*20+0.5,c(seq(1,length(wave_start_ix))),col="red")
text(as.Date(wave_dates_nondecimal$reset),growth_df$growth[wave_reset_ix]*20+0.5,c(seq(1,length(wave_reset_ix))),col="darkgreen")

#lines(cases_df$date,cases_df$cases/10000,col="green")

#growth_rate = TTR::ROC(m$fitted.values, type = "discrete")
#plot(date,growth_rate)
