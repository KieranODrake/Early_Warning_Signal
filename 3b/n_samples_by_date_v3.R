#' Exploring number of samples per day
#' Want to set the minimum cluster size to be a proportion of samples over 
#' past 8-12 weeks, so the parameter would vary over time with sequencing intensity

#' Version 3 uses the most recent sample date (after filtering for UK Pillar 2,  
#' removing erroneous samples an only including up to 99.9% quantile) as the end date.
#' Previously the date in the tree filename was used as the end date. 

#' Using the adjusted meta data to calculate the number of samples per day and over an 
#' N week period prior to the tree date
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02")
load("tfps-env-HPC-start_v5.RData")

#' Want the minimum descendants to be at least 20 regardless of the number of samples
#' So
#' Calculate time series of samples per day
library(data.table)
dt = as.data.table( amd_adj$sample_date )
sample_dt = dt[, .N , by = V1 ]
n_samples_df = data.frame(  "date" = sample_dt[ , 1 ]
                            , "n_samples" = sample_dt[ , 2 ])
colnames( n_samples_df ) = c( "date" , "n_samples" )
n_samples_df = n_samples_df[ order( n_samples_df[ , "date" ] ) , ] 
rm( dt )

par( mfrow = c( 3 , 1 ) )
plot(   n_samples_df$date
      , n_samples_df$n_samples
      , xlab = "date"
      , ylab = "number of samples"
      , main = "Number of SARS-CoV-2 samples")

#min(sample_dt[,2])


#' Obtain the dates for the most recent samples in each tree
 
check_df = data.frame( "tree_date" = replicate(length( tre_list ) , "NA" )
                       , "max_sample_freq" = replicate(length( tre_list ) , "NA" )
                       , "most_recent_sample" = replicate(length( tre_list ) , "NA" )
                       , "diff_days" = replicate(length( tre_list ) , "NA" )
                       , "max_date" = replicate(length( tre_list ) , "NA" )                         
)

sample_dates_list = list()

#' Read dates of phylogenetic trees
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/Trees/")
tree_files = list.files()
tre_dates <- as.Date( regmatches( tree_files , regexpr('(\\d{4}-\\d{2}-\\d{2})', tree_files) ) )

library(ape)
library(lubridate)

#' Process trees
for(i in 1:length(tre_list)){
  
  #' #' Read in tree
  #' Record names of files not able to be read
  tre_read_test <- suppressWarnings( try( tre <- ape::read.tree( tre_list[ i ] ) , silent = TRUE ) )
  if ( is.null( tre_read_test ) ) {    #class( tre_read_test ) == "try-error" ) {
    tre_read_error_list <- c( tre_read_error_list , tre_list[ i ] )
    message( "Error reading Tree")
    next
  } else { 
    tre = tre_read_test
    rm(tre_read_test) 
  }
  
  #' Trees need to be unrooted in order for the tfp scanner to work
  tre <- ape::unroot( tre )
  
  #' Check if tre$tip.label and amd$sequence_name format matches
  tre$tip.label         <- gsub( tre$tip.label         , patt = "\'", rep = "")
  amd_adj$sequence_name <- gsub( amd_adj$sequence_name , patt = "\'", rep = "")
  root_on_tip = "Wuhan/WH04/2020"
  
  #' First need to convert from tree data type to a dataframe data type
  #testdf = data.frame( stringr::str_sub( tre$tip.label , 2 , -2 ) ) # Tree tip labels are in '' whereas meta data is not
  testdf = data.frame( tre$tip.label )
  colnames(testdf) <- c("tre.tip.label")
  #' Define 'not in'
  '%ni%' <- Negate("%in%")
  #' Create lists of tip.labels that are in/not in the metadata
  exclude = subset( testdf , tre.tip.label %ni% amd_adj$sequence_name ) # not in the metadata (global tree vs UK metadata)
  include = subset( testdf , tre.tip.label %in% amd_adj$sequence_name ) # in the metadata (global tree vs UK metadata)
  if ( nrow( include ) == 0 ){
    tre_no_tips_in_meta_list <- c( tre_no_tips_in_meta_list , tre_list[ i ] )    
    next
  } 
  
  #' Remove non-UK tips from tree but keep the original 'Wuhan/WH04/2020'
  tre_UK_only = drop.tip( tre , exclude$tre.tip.label ) #[-c(1)]) # assumes that 'Wuhan/WH04/2020' is the first row in dataframe
  
  #' Remove NA sample dates
  tre_UK_only = drop.tip( tre_UK_only , no_sample_date$sequence_name )
  
  #' Run tfpscanner (emvolz-phylodynamics-tfpscan-v2_updated.R) for rolling periods
  #' and collect time series growth statistics
  
  #' Select 'root_on_tip' from Tree to be scanned
  tips <- subset( amd_adj, amd_adj$sequence_name %in% tre_UK_only$tip.label)
  root_on_tip <- subset( tips , tips$sample_date == min(tips$sample_date) )$sequence_name
  root_on_tip_sample_time <- decimal_date( min( tips$sample_date ) )
  
  # Remove unnecessary objects
  if ( exists( "test" ) ) { rm(test) }
  if ( exists( "tre" ) ) { rm(tre) }
  if ( exists( "include" ) ) { rm(include) }
  if ( exists( "exclude" ) ) { rm(exclude) }
  if ( exists( "amd" ) ) { rm(amd) }
  if ( exists( "tfps_output" ) ) { rm( tfps_output ) }
  if ( exists( "test_df" ) ) { rm( test_df ) }
  if ( exists( "tips" ) ) { rm( tips ) }
  gc()
  
  
  # Check frequency of sample names and return maximum frequency
  tips_table = table(data.frame(tre_UK_only$tip.label))
  #View(tips_table)
  check_df[ i, 2 ] = max( tips_table )
  
  #' Add tree date to results
  #' Create list of sample dates for each tree date
  sample_dates = subset( amd_adj , sequence_name %in% tre_UK_only$tip.label )$sample_date
  sample_dates_list[[i]] = sample_dates
  names( sample_dates_list[[ i ]] ) = tre_dates[i]
  
  #check_df[ i, 1 ] = as.Date( lubridate::ymd( data.frame( strsplit( tre_list[ i ] , "_" ) )[ 3 , 1 ]) , origin="1970-01-01" ) #extract date from filename
  check_df[ i, 1 ] = as.Date( data.frame( strsplit( tre_list[ i ] , "_" ) )[ 3 , 1 ] , origin="1970-01-01" ) #extract date from filename
  
  #' Find date of most recent sample in tree
  most_recent_date = max( subset( amd_adj , sequence_name %in% tre_UK_only$tip.label )$sample_date )
  
  check_df[ i, 3 ] = as.Date( most_recent_date , origin = "1970-01-01" )
  
  #' Compute difference between the tree date and the date of the most recent sample
  check_df[ i, 4 ] = as.numeric( check_df[ i, 1 ] ) - as.numeric( check_df[ i, 3 ] )
  
  message("Tree ",i," complete. Max freq: ",check_df[ i, 2 ],". Tree date: ",check_df[ i, 1 ],". Most recent sample date: ",check_df[ i, 3 ],". Date difference: ", check_df[ i, 4 ] )
}

n_samples_period_56_list = c()
n_samples_period_84_list = c()
#' Set period over which to sum the number of samples
for (i in 1 : length( tre_dates ) ){
  tre_end = as.Date( as.numeric( check_df[ i, 3 ] ) ,origin="1970-01-01") #as.Date( tre_dates[ i ] ) #"2022-03-29" )
  tre_start_56 = tre_end - 56 #' days / 8wks
  tre_start_84 = tre_end - 84 #' days / 12wks
  #' Calculate number of samples between tre_start and tre_end dates
  n_samples_period_56 = sum( subset( n_samples_df , ( n_samples_df$date >= tre_start_56 ) & ( n_samples_df$date <= tre_end ) )[ , 2 ] )
  n_samples_period_84 = sum( subset( n_samples_df , ( n_samples_df$date >= tre_start_84 ) & ( n_samples_df$date <= tre_end ) )[ , 2 ] )
  
  n_samples_period_56_list = c( n_samples_period_56_list , n_samples_period_56 )
  n_samples_period_84_list = c( n_samples_period_84_list , n_samples_period_84 )
}

n_samples_period_df = data.frame(   "date" = tre_dates
                 , "n_samples_period_56" = n_samples_period_56_list
                 , "n_samples_period_84" = n_samples_period_84_list
                 )

#' Plot sum of samples for two periods
plot(     n_samples_period_df$date
        , n_samples_period_df$n_samples_period_84
        , xlab = "date"
        , ylab = "number of samples"
        , xlim = c( min( n_samples_df$date ) , max( n_samples_df$date ) )
        , typ="l"
        , main = "Sum of samples over running 8 and 12 week periods"
        )
lines(  n_samples_period_df$date
      , n_samples_period_df$n_samples_period_56
      , col="red"
      )
legend(   "topleft"
        , legend = c( "Sum over 12 wks" , "Sum over 8 wks" )
        , col = c( "black" , "red" ) , lty = c( 1 , 1 ) )

#' Calculate the percentage in order to have 
#' the minimum number of cluster descendants as equal to 20
min_desc_perc_56 = 20 / min( n_samples_period_df$n_samples_period_56[-81] ) # row 81 removed as this relates to a tree that can't be read. If this is not done then the result is x/0=Inf
min_desc_perc_84 = 20 / min( n_samples_period_df$n_samples_period_84[-81] ) # row 81 removed as this relates to a tree that can't be read. If this is not done then the result is x/0=Inf

#' Create data frame of the min descendants calculated as a proportion of the
#' sum of samples over a period of preceding N days
min_desc_prop_df = data.frame(   "date" = n_samples_period_df$date
                               , "period_56" = round( n_samples_period_df$n_samples_period_56 * min_desc_perc_56 )
                               , "period_84" = round( n_samples_period_df$n_samples_period_84 * min_desc_perc_84 )
                               )
#' Save data frame of minimum descendants calculated as a proportion of the samples over previous 56 and 84 days
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/Trees/")
saveRDS( min_desc_prop_df , "min_desc_prop_df_v3.RData" )
#' Plot minimum descendants to be used in TFP Scanner for two periods
plot(     n_samples_period_df$date
          , n_samples_period_df$n_samples_period_84 * min_desc_perc_84
          , xlab = "date"
          , ylab = "minimum number of cluster descendants"
          , xlim = c( min( n_samples_df$date ) , max( n_samples_df$date ) )
          #, ylim = c( 0, 2000) 
          , ylim = c( 0 , max( max(n_samples_period_df$n_samples_period_84 * min_desc_perc_84) , max(n_samples_period_df$n_samples_period_56 * min_desc_perc_56) ))
          , typ="l"
          , main = "Min number of descendants for TFP Scanner calculated as % of sum of samples (based on absolute min of 20 descendants)"
          , col="blue"
          )
lines(  n_samples_period_df$date
        , n_samples_period_df$n_samples_period_56 * min_desc_perc_56
        , col="red"
)
legend(   "topleft"
          , legend = c( "12 wks" , "8 wks" )
          , col = c( "black" , "red" ) , lty = c( 1 , 1 ) )
max( n_samples_period_df$n_samples_period_56 * min_desc_perc_56 )

#' Plot cluster size for previous scans on top of minimum descendants to be used 
#' in TFP Scanner for two periods 
plot(     n_samples_period_df$date
          , n_samples_period_df$n_samples_period_84 * min_desc_perc_84
          , xlab = "date"
          , ylab = "minimum number of cluster descendants"
          , xlim = c( min( n_samples_df$date ) , max( n_samples_df$date ) )
          , ylim=c(0,2000)
          #, ylim = c( 0 , max( max(n_samples_period_df$n_samples_period_84 * min_desc_perc_84) , max(n_samples_period_df$n_samples_period_56 * min_desc_perc_56) ))
          , typ="l"
          , main = "Min descendants calculated as % of sum of samples (lower bound of 20 descendants)"
          , cex = 2
          )
lines(  n_samples_period_df$date
        , n_samples_period_df$n_samples_period_56 * min_desc_perc_56
        , col="red"
)

#' Add all the individual cluster sizes by date (multiple values for each date)
n = 18 #' number of the variable list eg 1 = min age 7, max age 56, min desc 20
for (j in 1 : 288 ){
  if ( n_samples_period_df$date[ j ] == "2021-01-01" ){ next }
  list_df = data.frame( tfps_cluster_size_lists[[n]][ j ] )
  if( nrow( list_df ) == 0 ){ next }
  for (k in 1 : nrow( list_df ) ){
    points(  subset( n_samples_period_df$date[ j ] , n_samples_period_df$date[ j ] != "2021-01-01" )
           , list_df[ k, ]
          , pch = 3
          , col = "grey"
          )
  }
}

lines(  n_samples_period_df$date
        , n_samples_period_df$n_samples_period_56 * min_desc_perc_56
        , col="red"
)

lines(  n_samples_period_df$date
        , n_samples_period_df$n_samples_period_84 * min_desc_perc_84
        , col="black"
)

legend(   "topleft"
          , legend = c( "12 wks" , "8 wks" , "cluster sizes for min age = 28, max age = 84, min desc = 100")
          , col = c( "black" , "red" , "grey" ) , lty = c( 1 , 1 ,NA ), pch=c(NA,NA,3) 
          #, cex = 2
          )


##########################################
#' Condensed version
#' Calculate minimum descendants as a proportion of samples in period (equal to maximum cluster age)
dt = as.data.table( amd_adj$sample_date )
sample_dt = dt[, .N , by = V1 ]
n_samples_df = data.frame(  "date" = sample_dt[ , 1 ]
                            , "n_samples" = sample_dt[ , 2 ])
colnames( n_samples_df ) = c( "date" , "n_samples" )
n_samples_df = n_samples_df[ order( n_samples_df[ , "date" ] ) , ] 
if ( exists( "dt" ) ){ rm( dt ) }
if ( exists( "sample_dt" ) ){ rm( sample_dt ) }

tre_end = check_df[ , 3 ] #as.Date( regmatches( tre_list[ as.numeric( params_run ) ] , regexpr('(\\d{4}-\\d{2}-\\d{2})', tre_list[ as.numeric( params_run ) ] ) ) )
tre_start = tre_end - period_length

n_samples_period = sum( subset( n_samples_df , ( n_samples_df$date >= tre_start ) & ( n_samples_df$date <= tre_end ) )[ , 2 ] )

if ( period_length == lubridate::days(56) ){ min_descendants = round( 0.003222688 * n_samples_period ) } #' multiplication factor calculated separately in 'n_samples_by_date.R' 
if ( period_length == lubridate::days(84) ){ min_descendants = round( 0.002373887 * n_samples_period ) } #' multiplication factor calculated separately in 'n_samples_by_date.R'
