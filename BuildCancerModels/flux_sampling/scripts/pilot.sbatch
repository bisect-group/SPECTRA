#!/bin/bash
#SBATCH --job-name=rev2_pilot
#SBATCH --account=def-mahadeva
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:30:00
#SBATCH --output=/home/pavan/scratch/SPECTRA_review_2/logs/pilot_%j.out
#SBATCH --error=/home/pavan/scratch/SPECTRA_review_2/logs/pilot_%j.err

module purge
module load matlab/2024b.1

export MODEL_ID=ACH-000681
export CHAIN=1
export NSAMPLES=500
export MAXTIME=2400
export OUTDIR=/home/pavan/scratch/SPECTRA_review_2/sampling/pilot

cd /home/pavan/scratch/SPECTRA_review_2
matlab -nodisplay -nosplash -nodesktop -r "\
  try, run('/home/pavan/scratch/SPECTRA_review_2/scripts/s05_sample_chain.m'); \
  catch e, disp(getReport(e)); exit(1); end; exit(0);"
