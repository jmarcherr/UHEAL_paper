#!/bin/sh
#BSUB -J run_preproc_abr
#BSUB -q hpc
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o log/Output_ABR.txt
#BSUB -e log/Error_ABR.txt
#BSUB -W 24:00 
#BSUB -n 1
#BSUB -u jonmarc@dtu.dk
#BSUB -N "Job done"
# -- end of LSF options --

cd //work3/jonmarc/UHEAL_paper/_eeg/_ABR/_private
matlab -nodisplay -r run_ABR_preproc_trials -logfile log/log_preproc_abr

