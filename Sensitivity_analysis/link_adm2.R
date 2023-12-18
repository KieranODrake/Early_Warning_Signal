#' Combining COG-UK metadata in order to match tree metadata to adm2 (region) data
#' Author: Kieran Drake
#' Date: October 2023

library(data.table)
library(dplyr)
library(sqldf)
library(stringr)

#' **Run first part of 1_tree_prep_pre_tfpscan_v5_adm1_Tg6_5.R to get the metadata used for the phylo trees ('amd_adj')**
#rm( bad_dates, no_sample_date, samples_to_be_removed_unique_name )
#' Add central sample id to amd_adj
split_seq_name = strsplit( amd_adj$sequence_name , "/" ) 
central_sample_id = sapply( split_seq_name ,"[", 2 )
amd_adj$central_sample_id = central_sample_id
rm( split_seq_name , central_sample_id )

#' Load meta data files containing amd2 data
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/linking_adm2/COG_adm2_Kieran/COG_adm2_Kieran/majora.latest.metadata.tsv/")
majora_latest_metadata <- fread("majora.latest.metadata.tsv")

#' Trim metadata to required columns
majora_latest_metadata_trimmed <- majora_latest_metadata[,c(6,2,3,4,8)]
majora_latest_metadata_trimmed$collection_date <- as.Date( majora_latest_metadata_trimmed$collection_date )
#' Check data
table(majora_latest_metadata$adm0)
table(majora_latest_metadata$adm1)
table(majora_latest_metadata$adm2)

#' Trim majora_latest_metadata_trimmed to sequences in amd_adj
#' match on central_sample_id
majora_latest_metadata_trimmed_2 = subset( majora_latest_metadata_trimmed 
                                         , majora_latest_metadata_trimmed$central_sample_id %in% amd_adj$central_sample_id )
#' Where more than one entry for a particular central_sample_id need to remove/merge additional rows, whilst retaining all information
table(duplicated( majora_latest_metadata_trimmed_2$central_sample_id ) )
table( majora_latest_metadata_trimmed_2$central_sample_id )
View( majora_latest_metadata_trimmed_2[ duplicated( majora_latest_metadata_trimmed_2$central_sample_id ),] )
majora_latest_metadata_trimmed_3 <- majora_latest_metadata_trimmed_2 %>% 
                                      dplyr::group_by( central_sample_id , adm1 , adm2 , adm2_private , collection_date ) %>% 
                                        summarise( na.rm = TRUE )
majora_latest_metadata_trimmed_3$na.rm = NULL
#' Check data frames
'%ni%' = Negate('%in%')
View( subset( amd_adj, amd_adj$central_sample_id %ni% majora_latest_metadata_trimmed_3$central_sample_id ) ) #Check - only central sample id not in majora metadata is the Wuhan sample that was manually entered into tree metadata.
View( subset( majora_latest_metadata_trimmed_3, majora_latest_metadata_trimmed_3$central_sample_id =="" ) )
View( subset( majora_latest_metadata_trimmed_3, majora_latest_metadata_trimmed_3$adm1 =="" ) )
View( subset( majora_latest_metadata_trimmed_3, majora_latest_metadata_trimmed_3$adm2 =="" ) ) #' 34,117 rows
View( subset( majora_latest_metadata_trimmed_3, majora_latest_metadata_trimmed_3$adm2_private =="" ) ) #' 351,110 rows
View( subset( majora_latest_metadata_trimmed_3, majora_latest_metadata_trimmed_3$adm2 =="" & majora_latest_metadata_trimmed_3$adm2_private =="" ) ) #' 18,997 rows
#'** There are 18,997 rows out of 2,155,743 where there is no adm2 or adm2_private (first part of postcode)**

#' Split 'majora_latest_metadata_trimmed_3' into rows with a value for adm2 or adm2_private and rows with no value for both
majora_latest_metadata_trimmed_3_val_no = subset( majora_latest_metadata_trimmed_3
                                                  , majora_latest_metadata_trimmed_3$adm2 =="" & majora_latest_metadata_trimmed_3$adm2_private =="" ) #' 18,997 rows
#adm2_meta_1_merge_val = subset(adm2_meta_1_merge , adm2_meta_1_merge %!in% adm2_meta_1_merge_no_val )
majora_latest_metadata_trimmed_3_val = subset( majora_latest_metadata_trimmed_3
                                               , majora_latest_metadata_trimmed_3$adm2 !="" | majora_latest_metadata_trimmed_3$adm2_private !="" ) #' 2,136,746 rows
#' Results in 18,997 (no values) + 2,136,746 (at least 1 value out of 2 categories) = 2,155,746, which matches the total in the original df

#' Rename so easier in later manipulation
majora_meta_val_1 = majora_latest_metadata_trimmed_3_val
majora_meta_val_no_1 = majora_latest_metadata_trimmed_3_val_no

#' Tidy up environment
rm( majora_latest_metadata , majora_latest_metadata_trimmed , majora_latest_metadata_trimmed_2 )
gc()

#' Need to clean adm2 geo data as there are some errors: 
#' - different names for same adm2 e.g. 'WINDSOR AND MAIDENHEAD' vs 'WINDSOR_AND_MAIDENHEAD' 
#' - London included as well as individual boroughs. Similar for other regions.
#' - some adm2 missing, so use postcode (if available to ascertain amd2)
#' Link to instructions from COG-UK on geo cleaning in UK https://github.com/COG-UK/docs/blob/master/geography_cleaning.md
#' Link to files for geo cleaning https://github.com/COG-UK/geography_cleaning/tree/master/geography_utils
#' First read in cleaning files and merge relevant columns with majora metadata
#' Second clean adm2 in majora metadata
#' Third add clean(est) adm2 to amd_adj

#' Read in cleaning files
setwd("C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_10/linking_adm2/")

#' Correct outer postcode mistakes
outer_postcode_cleaning <- fread("outer_postcode_cleaning.csv")
majora_meta_val_2_postcode <- merge( x = majora_meta_val_1 , y = outer_postcode_cleaning
                                   , by.x = "adm2_private" , by.y = "input_postcode"
                                   , all.x = TRUE 
                                   , all.y = FALSE )

#' Add column which uses the amended postcode or uses the Majora postcode if there is no amended postcode. Check using 'AG3' > 'SG3'
majora_meta_val_2_postcode$true_postcode_adj <- ifelse( is.na( majora_meta_val_2_postcode$true_postcode )
                                                        , majora_meta_val_2_postcode$adm2_private
                                                        , majora_meta_val_2_postcode$true_postcode)

#' Find adm2 from postcode
postcode_to_adm2 <- fread("postcode_to_adm2.tsv")
majora_meta_val_3_pc_to_adm2 <- merge( x = majora_meta_val_2_postcode , y = postcode_to_adm2
                                      , by.x = "true_postcode_adj" , by.y = "outer_postcode"
                                      , all.x = TRUE 
                                      , all.y = FALSE )
#' **3 rows added during the merge above**

#' We now have two adm2 columns so need to change their names
old_name_1 = "adm2.x" ; new_name_1 = "majora_adm2"
names( majora_meta_val_3_pc_to_adm2 )[ names( majora_meta_val_3_pc_to_adm2 ) == old_name_1 ] <- new_name_1
old_name_2 = "adm2.y" ; new_name_2 = "adm2_from_postcode"
names( majora_meta_val_3_pc_to_adm2 )[ names( majora_meta_val_3_pc_to_adm2 ) == old_name_2 ] <- new_name_2

#' Merge on LAD with adm2 in case LAD entered mistakenly instead of adm2
LAD_UTLA_adm2 <- fread("LAD_UTLA_adm2.csv") #' There are some rows with no value for LAD_name so need to remove them before merge
LAD_UTLA_adm2_trim_LAD <- subset( LAD_UTLA_adm2 , LAD_UTLA_adm2$LAD_name != "" )[,c("LAD_name","adm2","aggregated_adm2")]
majora_meta_val_4_LAD <- merge( x = majora_meta_val_3_pc_to_adm2 , y = LAD_UTLA_adm2_trim_LAD
                              , by.x = "majora_adm2" , by.y = "LAD_name"
                              , all.x = TRUE 
                              , all.y = FALSE )
#' **5 rows added during the merge above **

#' Change column names as going to make a second merge with 'LAD_UTLA_adm2' dataframe
names( majora_meta_val_4_LAD )[ names( majora_meta_val_4_LAD ) == "adm2"            ] <- "LAD_merge_adm2"
names( majora_meta_val_4_LAD )[ names( majora_meta_val_4_LAD ) == "aggregated_adm2" ] <- "LAD_merge_aggregated_adm2"

#' Once happy with meta dataframe, remove intermediate meta dataframes to save memory
rm( majora_meta_val_1 , majora_meta_val_2_postcode ,  majora_meta_val_3_pc_to_adm2 ) #, adm2_meta_5_UTLA )
gc()
#'** **

#' Change order of columns
majora_meta_val_4_LAD = majora_meta_val_4_LAD[ , c("central_sample_id"
                                                   ,"collection_date"
                                                   ,"adm1"
                                                   ,"majora_adm2"
                                                   ,"adm2_private"
                                                   ,"true_postcode"
                                                   ,"true_postcode_adj"
                                                   ,"adm2_from_postcode"
                                                   ,"LAD_merge_adm2"
                                                   ,"LAD_merge_aggregated_adm2")]

#' Replace incorrect (or inconsistent format) adm2 with correct version
adm2_cleaning <- fread("adm2_cleaning_w_col_heads.csv") #same as 'adm2_cleaning.tsv' file but with column headings added to all columns so can be easily read into R
#' There is a duplicate row for 'West Midlands' = 'WESTMIDLANDS' (no space) and 'GUERNSEY_CHANN' also duplicated, so needs to be remove both
length( unique( adm2_cleaning$location_provided_in_input ) ) #'199 unique values and 200 rows
View( adm2_cleaning[ duplicated( adm2_cleaning ) , ] )
adm2_cleaning_trim = subset( adm2_cleaning , adm2_cleaning$location_provided_in_input != "WESTMIDLANDS" )
adm2_cleaning_trim[which(adm2_cleaning_trim$location_provided_in_input == "GUERNSEY_CHANNE")]
adm2_cleaning_remove = subset( adm2_cleaning_trim , adm2_cleaning_trim$location_provided_in_input == "GUERNSEY_CHANNE" )
adm2_cleaning_keep = subset( adm2_cleaning_trim , adm2_cleaning_trim$location_provided_in_input != "GUERNSEY_CHANNE" )
adm2_cleaning_trim_2 = rbind( adm2_cleaning_keep , adm2_cleaning_remove[1,])
#adm2_cleaning_trim[199,] = NULL

majora_meta_val_5_ac <- merge( x = majora_meta_val_4_LAD , y = adm2_cleaning_trim_2
                              , by.x = "majora_adm2" , by.y = "location_provided_in_input"
                              , all.x = TRUE 
                              , all.y = FALSE )

#' Create definitive list of adm2 regions
adm2_complete_list = unique(LAD_UTLA_adm2_trim_LAD$adm2)

#' See 'adm2 decision tree.ppt' for visual representation of the following code

#' If the adm2 linked to the postcode is in the complete list of adm2s 
#' (i.e. it has a value and it multiple adm2s are not listed)
#' then set as 'adm2_final' (because it is determined from the postcode which is 
#' more likely to be correct (after cleaning))
majora_meta_val_6 = majora_meta_val_5_ac
majora_meta_val_6$adm2_final = apply( majora_meta_val_6 , 1 , function(x){
                                      if ( x["adm2_from_postcode"] %in% adm2_complete_list ){ #' Check if it is a real adm2, i.e. not NA, not multiple adm2s, and not an erroneous adm2
                                        result = x["adm2_from_postcode"] #  x["adm2_final"] == x["adm2_from_postcode"] 
                                      } else { result = "" }
                                      return( result )
                                    })

#' Next option for adm2_final' is the 'LAD_merge_adm2' which corrects the 'majora_adm2' from a LAD to an adm2 if a LAD was entered by mistake
majora_meta_val_6$adm2_final = apply( majora_meta_val_6 , 1 , function(x){
                                    if( x["adm2_final"] == "" ){
                                      if ( x["LAD_merge_adm2"] %in% adm2_complete_list ){ #' Check if it is a real adm2, i.e. not NA, not multiple adm2s, and not an erroneous adm2
                                        result = x["LAD_merge_adm2"]
                                      } else { result = x["adm2_final"] } #' If already a value for adm2_final then don't change it
                                    } else { result = x["adm2_final"] } #' If already a value for adm2_final then don't change it
                                    return( result )
                                  })

#' Next option for adm2_final' is the 'correct_location_1' if there are not multiple
#' "correct_locations". This is a correction for incorrect 'majora_adm2' values i.e. misspellings etc
majora_meta_val_6$adm2_final = apply( majora_meta_val_6 , 1 
                                         , function(x){
                                                        if( x["adm2_final"] == "" ) {
                                                          if( is.na(x["correct_location_2"]) | x["correct_location_2"] =="" ){ #' Make there are not multiple correct locations
                                                            if ( x["correct_location_1"] %in% adm2_complete_list ){ #' Check if it is a real adm2, i.e. not NA, and not an erroneous adm2
                                                              result = x["correct_location_1"] 
                                                            } else { result = x["adm2_final"] } #' If already a value for adm2_final then don't change it
                                                          } else { result = x["adm2_final"] } #' If already a value for adm2_final then don't change it
                                                        } else { result = x["adm2_final"] } #' If already a value for adm2_final then don't change it
                                                        return( result )
                                                      }
                                        )
#' Next option is to assume that 'majora_adm2' is correct (as long as it is a real adm2 i.e. in the 'adm2_complete_list')
majora_meta_val_6$adm2_final = apply( majora_meta_val_6 , 1 
                                         , function(x){
                                           if( x["adm2_final"] == "" ) {
                                             if( !is.na(x["majora_adm2"]) & x["majora_adm2"] !="" ){ #' Make sure there is a value
                                               if ( x["majora_adm2"] %in% adm2_complete_list ){ #' Check if it is a real adm2, i.e. not NA, and not an erroneous adm2
                                                 result = x["majora_adm2"]
                                               } else { result = x["adm2_final"] } #' If already a value for adm2_final then don't change it
                                             } else { result = x["adm2_final"] } #' If already a value for adm2_final then don't change it
                                           } else { result = x["adm2_final"] } #' If already a value for adm2_final then don't change it
                                           return( result )
                                         }
)

#' 279,259 rows still without adm2 plus the 18,997 rows split out earlier
nrow( subset( majora_meta_val_6 , majora_meta_val_6$adm2_final == "" ) )
majora_meta_clean = subset( majora_meta_val_6 , majora_meta_val_6$adm2_final !="" )
nrow( majora_meta_clean ) #' 1,857,495 sequences with clean adm2 (1,857,487 = 2,155,743 - 279,259 - 18,997)
majora_meta_clean = majora_meta_clean[,c("central_sample_id","collection_date","adm1","adm2_final")]
View(table( majora_meta_clean$adm2_final ))

#' majora_meta_clean can then be used in '1_tree_prep_pre_tfpscan_v5_adm2.R'

#' Further work could be done to try to assign adm2 values to the remaining sequences