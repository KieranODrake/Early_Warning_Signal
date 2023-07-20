#' Takes data frames created on HPC using top_variants_lgr_xxx_pval_xxx.R. These dataframes represent the 
#' top five variants by cluster frequency and growth rate for each phylogenetic tree and each set of cluster
#' parameters:
#'  - minimum cluster size
#'  - maximum cluster size
#'  - minimum descendants
#' and filter variables:
#'  - logistic growth rate (LGR) p-value
#'  - parent/sub-cluster replacement LGR level
#'   
#' This code takes all of the .rds files in a folder and combines them into a single dataframe, which can 
#' then be referenced when assessing early warning signals (EWS) as either True or False Positives.

#' Read in dataframe files, merge dataframes and save merged dataframe
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/top_variants/dataframes_output/" )
file_list = list.files() ; file_list = subset( file_list , file_list != "desktop.ini" )
top_variants_list = lapply( file_list , readRDS )
top_variants_df = data.table::rbindlist( top_variants_list )
setwd( "C:/Users/kdrake/OneDrive - Imperial College London/Documents/TFPS/tfps runs/2023_02/quantile_historical/top_variants/" )
saveRDS( top_variants_df , "top_variants_df.rds")
top_variants_df =readRDS("top_variants_df.rds")
#' Check that all TFP Scan parameters and filter variables have been included
#' Build list of parameters present in files generated on HPC
top_var_parameter_list = as.data.frame( paste0(   top_variants_df$TFPS_cluster_min_age , "_" 
                                                , top_variants_df$TFPS_cluster_max_age , "_" 
                                                , top_variants_df$TFPS_cluster_min_descendants , "_" 
                                                , top_variants_df$LGR_p_value_threshold , "_" 
                                                , top_variants_df$parent_sub_LGR_threshold 
                                               )
                                      )
top_var_parameter_list = unique( top_var_parameter_list )

#' Build list of expected parameters
min_age = c("07","14","28")
max_age = c("56","84")
min_desc = c("100","20","50","proportional")
p_val_th = c(0.01,0.05,10000)
parent_sub_lgr = c(0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1, NA)

check_list = data.frame()
for ( mina in min_age ){
  for ( maxa in max_age ){
    for ( mind in min_desc ){
      for ( pval in p_val_th ){
        for ( psub in parent_sub_lgr ){
          new_row = as.data.frame( paste0(   mina , "_" 
                                             ,  maxa , "_" 
                                             ,  mind , "_" 
                                             ,  pval , "_" 
                                             ,  psub
                                            )
                                    )
          check_list = rbind( check_list , new_row )
          
        }
      }
    }
  }
}

'%!in%' = Negate('%in%')
#nrow(subset( check_list , top_var_parameter_list %in% check_list ))
#top_var_parameter_list[1,1] == check_list[1,1]
#temp = check_list[ top_var_parameter_list[,1] %!in% check_list[,1] , ]
#setdiff(sort(check_list) , sort(top_var_parameter_list))
#setdiff(sort(top_var_parameter_list) , sort(check_list) )

library(sqldf)
sqldf('SELECT * FROM check_list EXCEPT SELECT * FROM top_var_parameter_list')

#' 57 of the expected 720 files (sets of parameters and filters) were not created as 
#' the HPC job ran out of memory (8gb wall set). These sets of parameters/filters 
#' were re-run with 16gb limit on HPC and worked.