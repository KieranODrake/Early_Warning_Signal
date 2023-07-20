#' Loads data into R for use in tfpscanner.R
#' Public data can be downloaded from: 
#' https://www.cogconsortium.uk/priority-areas/data-linkage-analysis/public-data-analysis/
#' Latest phylogenetic Tree available at 
#' https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_global_tree.newick
#' Using pre-downloaded phylogenetic trees in 
#' C:\Users\kdrake\OneDrive - Imperial College London\Documents\Transmission Fitness Polymorphism scanner (tfpscanner)\Trees
#' Metadata https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz

########################
# Load libraries
library("ape")
library("data.table")
library("R.utils")
#library("ggtree")
library("BiocManager")
library("lubridate")
#######################

#' To enable HPC array number from command line to be used as an input in this R script
params_run <- commandArgs(trailingOnly = TRUE)

#' Split when array parameter input file has multiple inputs
#teste2 <- str_split(params_run, pattern = ",")[[1]][1]

#'** Set parameter variables**
#' Minimum cluster age
min_cluster_age_yrs = 7/365 # in years # 14/365 # 28/365
#' Maximum cluster age
period_length = lubridate::days(56) # days (56 days = 8 weeks = 2 months = minimum for reasonable growth stats). Alternative is 84 days = 12 weeks
#' Minimum cluster size
min_descendants = 20 # 50 # 100


#' Load R environment saved with all the required functions (tfpscan), 
#' data (amd_adj and no_sample_date) and values (scan_error_list, tre_list,
#' tre_no_tips_in_meta_list, tre_read_error_list)
#setwd("Z:/home/tfp_scanner")
#setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/tfps runs/2022_09")
load("tfps-env-HPC-start_v5.RData")

#total_start_time = Sys.time()

#' No need for for loop if using HPC array jobs
#for ( m in 1 : length( params_run ) ){ #1 : 2 ){ # : length( tre_list ) ) {
#i = as.numeric( str_split( params_run , pattern = "," )[[ m ]][ 1 ] )
i = as.numeric( params_run ) #+ 2  #[[ m ]][ 1 ]

  scan_start_time = Sys.time()
  message(i," ", tre_list[ i ])
  
  #' Read in tree
  #' Record names of files not able to be read
  #setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Transmission Fitness Polymorphism scanner (tfpscanner)/Trees")
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

  # Set parameters for rolling time period
  #period_length = lubridate::days(84) # days (56 days is 8 weeks = 2 months = minimum for reasonable growth stats)
  phylo_tree_date = lubridate::ymd( data.frame( strsplit( tre_list[ i ] , "_" ) )[ 3 , 1 ]) #extract date from filename
  #end_date = subset( most_recent_sample_df , tree_file_date == phylo_tree_date )$most_recent_sample
  end_date = lubridate::ymd( max( tips$sample_date ) )
  start_date = end_date - period_length
  
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
  
  # Loop from start to finish of total time series
  #' Using emvolz-phylodynamics-tfpscanner.R 
  #' taken from Github https://github.com/emvolz-phylodynamics/tfpscanner/blob/master/R/tfpscanner.R last commit in Feb 2022

  tfps_output = tfpscan(  tre = tre_UK_only
                        , amd = amd_adj # = md #amd_adj
                        , min_descendants = min_descendants 
                        , max_descendants = 20e3
                        , min_cluster_age_yrs = min_cluster_age_yrs
                        , min_date = start_date # NULL
                        , max_date = end_date #NULL 
                        , min_blen = 1/30e3/2
                        , ncpu = 1
                        , output_dir = paste0('tfpscan-', Sys.Date())
                        , num_ancestor_comparison = 500
                        , factor_geo_comparison = 5
                        , Tg = 6.5/365
                        , report_freq = 50 
                        , mutation_cluster_frequency_threshold = 0.75
                        , test_cluster_odds = c() 
                        , test_cluster_odds_value = c() 
                        , root_on_tip = root_on_tip #"'Wuhan/WH04/2020'"
                        , root_on_tip_sample_time = root_on_tip_sample_time #2020 
                        , detailed_output = FALSE 
                        , compute_gam = TRUE
                        , compute_cluster_muts = FALSE
                        , tree_date = phylo_tree_date
                        #, compute_cluster_tree = FALSE
                        )

    scan_end_time = Sys.time()
    print( paste("scan ", i ,":") )
    print( paste( scan_end_time - scan_start_time ) )
    
#}
#total_end_time = Sys.time()
#print( paste( "Total time taken: ", (total_end_time - total_start_time) / 3600 , "hours" ) )
