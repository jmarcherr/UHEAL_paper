#!/bin/sh
#BSUB -J run_clin_scraper
#BSUB -q hpc
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o log/Output_clin.txt
#BSUB -e log/Error_clin.txt
#BSUB -W 24:00 
#BSUB -n 1
#BSUB -u jonmarc@dtu.dk
#BSUB -N "Job done"
# -- end of LSF options --

cd //zhome/7e/f/64621/Desktop/UHEAL_paper/_clin
matlab -nodisplay -r run_clin_scraper -logfile log/log_clin

