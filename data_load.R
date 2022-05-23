#' Function to load data for Covid-19 cases or hospitalisations
#' and to reformat
#' @param folder Location of data file.
#' @param filename Name of file to beloaded
#' @param dat_type Type of data: Covid-19 "cases" or "hospitalisations" 
#' @return Dataframe with reformatted data

data_load <- function( folder, filename, dat_type ) 
{
  library( lubridate ) 
  
  # Load data, trim and organise
  setwd(folder)
  cases_df = fread( filename )

  if (dat_type == "cases"){
        cases_df = cases_df[,c("date","newCasesBySpecimenDate")]
  }
  if (dat_type == "hospitalisations"){
        cases_df = cases_df[,c("date","newAdmissions")]
  }

  colnames(cases_df) <- c("date","cases") # Rename columns
  cases_df = cases_df[!is.na(cases_df$cases),]# Remove days with NA cases
  # Change order of data so it is ascending and organise dates
  # (otherwise the GAM derivative series is in reverse time order)
  if(cases_df$date[1]>cases_df$date[2]){
    cases_df$date <- rev(cases_df$date)
    cases_df$cases <- rev(cases_df$cases)
  }
  cases_df$date <- as.Date( cases_df$date ) # Change format of date
  cases_df$time <- lubridate::decimal_date(cases_df$date) # Add column for decimal date
  cases_df$wday <- lubridate::wday( cases_df$date ) # Add column for day of week

  par(mfrow=c(1,1))
  plot(cases_df$date,cases_df$cases,xlab = "Date", ylab = dat_type)
  return(cases_df) # Return from function
}