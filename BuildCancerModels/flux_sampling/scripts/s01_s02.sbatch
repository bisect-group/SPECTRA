#!/bin/bash
#SBATCH --job-name=rev2_build
#SBATCH --account=def-mahadeva
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/home/pavan/scratch/SPECTRA_review_2/logs/build_%j.out
#SBATCH --error=/home/pavan/scratch/SPECTRA_review_2/logs/build_%j.err

module purge
module load gurobi
module load matlab/2024b.1
export GUROBI_PATH="${GUROBI_HOME}"

cd /home/pavan/scratch/SPECTRA_review_2

matlab -nodisplay -nosplash -nodesktop -r "\
  try, \
    run('/home/pavan/scratch/SPECTRA_review_2/scripts/s01_prepare_depmap_data.m'); \
    clear; \
    run('/home/pavan/scratch/SPECTRA_review_2/scripts/s02_build_models.m'); \
  catch e, \
    disp(getReport(e)); exit(1); \
  end; exit(0);"
