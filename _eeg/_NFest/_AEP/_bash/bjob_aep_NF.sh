#!/bin/sh
#BSUB -J run_preproc_abr
#BSUB -q hpc
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o log/Output_AEP.txt
#BSUB -e log/Error_AEP.txt
#BSUB -W 24:00 
#BSUB -n 1
#BSUB -u jonmarc@dtu.dk
#BSUB -N "Job done"
# -- end of LSF options --

cd //work3/jonmarc/UHEAL_paper/_eeg/_NFest/_AEP
matlab -nodisplay -r run_AEP_preproc_NF -logfile log/log_preproc_aep

