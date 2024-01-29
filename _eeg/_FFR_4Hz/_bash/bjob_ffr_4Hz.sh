#!/bin/sh
#BSUB -J run_preproc_ffr_4hz
#BSUB -q hpc
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o log/Output_FFR_4Hz.txt
#BSUB -e log/Error_FFR_4Hz.txt
#BSUB -W 24:00 
#BSUB -n 1
#BSUB -u jonmarc@dtu.dk
#BSUB -N "Job done"
# -- end of LSF options --

cd //work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_scripts
matlab -nodisplay -r run_FFR_4Hz_preproc -logfile log/log_preproc_ffr_4Hz

