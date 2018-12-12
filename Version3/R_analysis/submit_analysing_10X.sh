#!/usr/bin/env bash 
 # @ cpus_per_task = 8 
 # @ job_name = 10X_R 
 # @ output = 10X_R.out 
 # @ error = 10X_R.err  
 
 module load R/3.5.0 


Rscript analysing_10Xdata.R
