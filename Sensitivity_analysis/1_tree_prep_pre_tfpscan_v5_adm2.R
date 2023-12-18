#' Process phylogenetic tree data prior to loading into tfpscanner.R
#' Public data can be downloaded from: 
#' https://www.cogconsortium.uk/priority-areas/data-linkage-analysis/public-data-analysis/
#' Latest phylogenetic Tree available at 
#' https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_global_tree.newick
#' Using pre-downloaded phylogenetic trees in 
#' C:\Users\kdrake\OneDrive - Imperial College London\Documents\Transmission Fitness Polymorphism scanner (tfpscanner)\Trees
#' Metadata https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz

########################
# Load packages and libraries
install.packages("data.table")
install.packages('R.utils')
install.packages ("ggplot")
install.packages("lubridate")
if (!require ("BiocManager", quietly = TRUE)) install.packages ("BiocManager") BiocManager::install ("ggtree")
install.packages("BiocManager") 
BiocManager::install ("ggtree")
library(ape)
library(data.table)
library(R.utils)
library(ggplot2)
library(ggtree)
library(lubridate)
#######################

#' Need to load tfpscan function 
#' emvolz-phylodynamics-tfpscanner.R taken from github with last commit in Feb 2022. 
#' KD amendments in rbind code in prune tree section to make columns the same for objects to be rbind and also output file names.

# Load sequence meta data to accompany phylogenetic trees
#meta_filename = file.choose()
#meta = fread(meta_filename)
#' Or
meta = fread("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/cog_metadata_2022_05_03.csv.gz")
#meta_csv= read.csv("cog_metadata_2022_05_03.csv.gz")

#' Put metadata into dataframe and only select required variables
amd = data.frame(   meta$sequence_name
                  , meta$sample_date
                  , meta$adm1
                  , meta$lineage
                  #, meta$mutations
                  , meta$is_pillar_2 
                  )
#amd$meta.mutations = "A" # set all mutations to "A" to speed up for run NOT CURRENTLY WORKING WITHOUT MUTATIONS

# Change amd dataframe column names to match tfpscan function 
colnames( amd ) <- c(   "sequence_name"
                      , "sample_date"
                      , "region"
                      , "lineage"
                      #, "mutations"
                      , "pillar_2")

# Add metadata for 'Wuhan/WH04/2020'
amd_adj = data.frame( amd )
amd_adj[ nrow( amd_adj ) + 1 , ] = c(   "Wuhan/WH04/2020"
                                      , "2020-01-05"
                                      , "CHINA"
                                      , "Original"
                                      #, "A" #"orf1ab:S135R") # Sample date from https://www.epicov.org/epi3/frontend#2155fc for EPI_ISL ID 406801
                                      , "Y")  #'Y so that is not removed when filter for Pillar 2

# Make sure data is in correct formats
amd_adj$sample_date = lubridate::ymd( amd_adj$sample_date )
amd_adj$sequence_name = as.character( amd_adj$sequence_name )
#amd_adj$country = as.character(amd_adj$country)
amd_adj$pillar_2 = as.character( amd_adj$pillar_2 )
amd_adj$region = as.character( amd_adj$region )
amd_adj$lineage = as.character( amd_adj$lineage )
#amd_adj$mutations = as.character( amd_adj$mutations )
#' Remove NA sample dates
no_sample_date = amd_adj[ is.na( amd_adj$sample_date ) , ]
amd_adj = amd_adj[ !is.na( amd_adj$sample_date ) , ]
#' Filter for Pillar 2 (community) samples
amd_adj = amd_adj[ amd_adj$pillar_2 == "Y" , ]
#' Then remove Pillar 2 column as no longer needed
amd_adj$pillar_2 = NULL
#' Remove sequences found to have sample dates after the relevant tree dates - see 'seq_tree_data_lag_calc.R' for code
bad_dates = data.frame(   "England/PHEC-228D3/2020"
                        , "England/PORT-2D9068/2021"
                        , "England/PORT-1E6E843/2021"
                        , "England/NORT-1BAD83C/2021"
                        , "England/PHEC-30AC96/2021"
                        , "England/NORW-3016752/2021"
                        , "England/PHEP-018655/2021"
                        , "England/PHEP-018200/2021"
                        , "England/PHEP-025619/2021"
                        , "England/PHEP-025615/2021"
                        , "England/NORW-303A602/2021"
                        , "England/NORW-303E7A4/2021"
                        , "England/NORW-3047E8F/2021"
                        , "England/NORW-304C1C9/2021"
                        , "England/NORW-306F974/2021"
                        , "England/NORW-306DC2F/2021"
                        )

'%ni%' <- Negate("%in%") #' Define 'not in'

amd_adj = amd_adj[ amd_adj$sequence_name %ni% bad_dates , ]

#' Remove sequences dated after the 99.9% quantile sample date - see 'seq_tree_data_date_error_calc.R' for code
#' Quantiles were calculated for each tree but samples falling in the 0.1% are to be removed for all trees
#' This adjustment is made to remove those trees that have dates unfeasibly close to the tree date i.e. there is
#' not enough time between sample date and tree date for genomic sequencing and phylogenetic tree pipeline
#' Note that the list below will include the 'bad_dates' list above
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/Early Warning Signal/Data lags/Possible sample date errors")
samples_to_be_removed_unique_name <- readRDS( "samples_to_be_removed_unique_name.rds" )
amd_adj = amd_adj[ amd_adj$sequence_name %ni% samples_to_be_removed_unique_name , ]
#' A total of 22,099 unique samples removed (including 16 with dates after the tree date) from a total of 2,177,827 samples.
#' This is equivalent to 1.01% of samples. The proportion removed from each 56 day or 84 day analysis period will be less.

# Check/Remove duplicates of patient, date and lineage (i.e. sequences relating to the same infection, but not repeat infections of same patient)
#identifier = paste0(amd_adj$sequence_name,"/",amd_adj$sample_date,"/",amd_adj$lineage)
#amd_adj$identifier = identifier
#freq_table = table(amd_adj$identifier)
#freq_table_order = freq_table[order(freq_table,decreasing=TRUE)]
#View(freq_table_order)
# No duplicates found in the metadata

# Remove to free up memory
rm( meta , amd )
gc()

#############################
#' ** Run through 'link_adm2.R' to find adm2 regions to link to the sequence names in 'amd_adj' **
#' 'link_adm2.R' output of interest is 'majora_meta_clean', which contains adm2 values
#' Merge amd_adj and majora_meta_clean on central_sample_id
amd_adj_merge <- merge( x = amd_adj , y = majora_meta_clean
                                     , by.x = "central_sample_id" , by.y = "central_sample_id"
                                     , all.x = TRUE 
                                     , all.y = FALSE )
nrow( subset( amd_adj_merge , is.na( amd_adj_merge$adm2_final ) ) ) / nrow( amd_adj_merge ) #' 13.8% without an adm2
#' Retain copy of old amd_adj before replacing with adm2 adj
amd_adj_adm1 = amd_adj 
amd_adj = amd_adj_merge[,c("sequence_name","sample_date","adm2_final","lineage")]
colnames(amd_adj) = c("sequence_name","sample_date","region","lineage")
amd_adj = subset( amd_adj , !is.na(amd_adj$region) )
View( table( amd_adj$region ) )
#' Check for duplicates

#' Add back in Wuhan meta data (removed in merger and adm2 NA filter)
amd_adj[ nrow( amd_adj ) + 1 , ] = c(   "Wuhan/WH04/2020"
                                        , "2020-01-05"
                                        , "CHINA"
                                        , "Original") 
                                        #, "A" #"orf1ab:S135R") # Sample date from https://www.epicov.org/epi3/frontend#2155fc for EPI_ISL ID 406801
amd_adj$region = as.character( amd_adj$region )
#' Look at timing of sequences without adm2

#' backup environment before removing variables ready for analysis
save.image( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/linking_adm2/manipulating_adm2_env.RData" )
#' Load environment as required
load("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/linking_adm2/manipulating_adm2_env.RData" )

rm(adm2_cleaning, adm2_cleaning_keep , adm2_cleaning_remove , adm2_cleaning_trim , adm2_cleaning_trim_2
   , amd_adj_merge , LAD_UTLA_adm2, LAD_UTLA_adm2_trim_LAD , majora_latest_metadata_trimmed_3
   , majora_latest_metadata_trimmed_3_val , majora_latest_metadata_trimmed_3_val_no , majora_meta_clean
   , majora_meta_val_1 , majora_meta_val_4_LAD , majora_meta_val_5_ac, majora_meta_val_6 , majora_meta_val_no_1
   , outer_postcode_cleaning , postcode_to_adm2 , adm2_complete_list , new_name_1 , new_name_2 , old_name_1
   , old_name_2 )
saveRDS( amd_adj_adm1 , "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/linking_adm2/amd_adj_adm1.rds" )
amd_adj_adm1 = readRDS( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/linking_adm2/amd_adj_adm1.rds" )

rm(amd_adj_adm1)
gc()

#############################

#' Select directory containing phylogenetic trees and create list
#tre_dir = choose.dir()
#setwd( tre_dir )
#' Or
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/Trees")

tre_list <- list.files()
tre_list #' Check they are all phylogenetic trees
tre_list <- tre_list[1:289] #' There are some other, non-tree, files saved in this directory now
tre_list #' Check they are all phylogenetic trees

# Initialise list of files not able to be read
tre_read_error_list <- c()
tre_no_tips_in_meta_list <- c()
scan_error_list <- c()

#' ** Load TFP Scanner **
#' C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/emvolz-phylodynamics-tfpscanner.R

#' Save environment for use in High Performance Computing (HPC) cluster
#save.image( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/tfps-env-HPC-start_v5.RData" )
save.image( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/tfps-env-HPC-start_v5_adm2.RData" )
#' Load environment as required
#load("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/tfps-env-HPC-start_v5.RData" )
load("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/tfps-env-HPC-start_v5_adm2.RData" )
