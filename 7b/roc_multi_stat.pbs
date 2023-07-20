#!/bin/bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=1:mem=8gb:gpfs=true
#PBS -J 1-30

## Load environment
module load anaconda3/personal
source ~/anaconda3/etc/profile.d/conda.sh
conda activate R4

## Create temporary directory on the computer node
export JOB_NUM=$(echo ${PBS_JOBID} | cut -f 1 -d '.' | cut -f 1 -d '[')
export WORKDIR="${EPHEMERAL}/${JOB_NUM}.${PBS_ARRAY_INDEX}"
mkdir -p $WORKDIR

## Copy across required files to job working directory
## First copy files from this directory so that have information for which scan directories to copy
cp -R $HOME/tfps_2023_02/quantile_historical/analysis/ROC_multi_stat/array_input_multi.txt $WORKDIR

## Move to Working directory
cd $WORKDIR

## Create variables for parent/sub-cluster LGR threshold level and LGR p-value for use in folder names
export PS_LGR=$( head -n ${PBS_ARRAY_INDEX} array_input_multi.txt | tail -1 | cut -d',' -f1)
export PVAL=$( head -n ${PBS_ARRAY_INDEX} array_input_multi.txt | tail -1 | cut -d',' -f2)

cp -R $HOME/tfps_2023_02/quantile_historical/analysis/ROC_multi_stat/roc_multi_stat.pbs $WORKDIR
cp -R $HOME/tfps_2023_02/quantile_historical/analysis/ROC_multi_stat/EWS_gen_time_window_TFPS_multi_stat_HPC.R $WORKDIR
cp -R $HOME/tfps_2023_02/quantile_historical/analysis/large_cluster_adjust_lgr_threshold_$PS_LGR/p_val_filter_$PVAL/dataframes_statistics/* $WORKDIR

## Move to Working directory
#cd $WORKDIR

## Run R script with job array number as input
Rscript $WORKDIR/EWS_gen_time_window_TFPS_multi_stat_HPC.R $(head -n $PBS_ARRAY_INDEX array_input_multi.txt | tail -1)

## Copy output files from job temporary working directory
cp -R $WORKDIR/ROC_multi_stat_*.rds $HOME/tfps_2023_02/quantile_historical/analysis/ROC_multi_stat/output/

conda deactivate
