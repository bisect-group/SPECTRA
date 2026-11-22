#!/bin/bash
#SBATCH --job-name=rev2_samp
#SBATCH --account=def-mahadeva
#SBATCH --array=1-80%40
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=5:00:00
#SBATCH --output=/home/pavan/scratch/SPECTRA_review_2/logs/samp_%A_%a.out
#SBATCH --error=/home/pavan/scratch/SPECTRA_review_2/logs/samp_%A_%a.err

module purge
module load matlab/2024b.1

R=/home/pavan/scratch/SPECTRA_review_2
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $R/scripts/sampling_tasks.txt)
export MODEL_ID=$(echo $LINE | cut -d' ' -f1)
export CHAIN=$(echo $LINE | cut -d' ' -f2)
export NSAMPLES=2500
export MAXTIME=16200          # 4.5h, inside the 5h wall limit
export OUTDIR=$R/sampling/chains

echo "task $SLURM_ARRAY_TASK_ID -> $MODEL_ID chain $CHAIN"
cd $R
matlab -nodisplay -nosplash -nodesktop -r "\
  try, run('$R/scripts/s05_sample_chain.m'); \
  catch e, disp(getReport(e)); exit(1); end; exit(0);"
