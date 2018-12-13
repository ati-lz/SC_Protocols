#!/usr/bin/env bash 
 # @ cpus_per_task = 8 
 # @ job_name = all_R 
 # @ output = all_R.out 
 # @ error = all_R.err  
 
 module load R/3.5.0 


Rscript analysing_alldata.R
