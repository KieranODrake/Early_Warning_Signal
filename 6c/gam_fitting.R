#' Function to fit generalised additive model to Covid-19 cases/hospitalisation data
#' Adapted by KD (May 2022) from sc2growth.R downloaded from 
#' https://gist.github.com/emvolz/3d858e03f3780b953067901121e3c7c0 on 4 May 2022
#' @param cases_df Dataframe containing case/hospitalisation data with columns: Date,cases,time,wday.
#' @param GAM_smooth_function Code for smoothing function in generalised 
#'    additive model from mgcv package. See help(smooth.terms) for descriptions.
#'    Options are: "tp","ts","ds","cr","cs","sos","ps","cp","re","gp" or "so". 
#' @param deg_free_k Degrees of freedom used for time element in GAM model fit.
#'    See help(choose.k) for further information.
#' @return m Generalised additive model (GAM)

gam_fitting <- function( cases_df, GAM_smooth_function = "cc", deg_free_k = 30 ) 
{
  library( lubridate )
  library( mgcv )
  
  # Fit generalised additive model (GAM)
  m <- mgcv::gam( cases ~ s( time, bs=GAM_smooth_function, k = deg_free_k) + s( wday, bs="cc", k = 7), data = cases_df ) #family = quasipoisson,
  #m <- mgcv::gam( cases ~ s( time, bs=GAM_smooth_function, k = deg_free_k ), family = quasipoisson, data = cases_df )
  
    # Plot graphs to check quality of model fit
  par(mfrow=c(2,2))
  #gam.check(m)
  #summary(m)
  return(m) # Return from function
}