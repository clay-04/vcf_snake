#!/bin/bash
#SBATCH --job-name=snakemake_bamvcf
#SBATCH --output=logs/snakemake_%j.log
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --partition=epyc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         # Just enough for Snakemake controller
#SBATCH --mem=12G                  # Enough for the controller itself
#SBATCH --time=240:00:00
#SBATCH --mail-user=clay007@ucr.edu
#SBATCH --mail-type=ALL

# Load environment
export PS1=""
source /etc/profile
module purge
module load snakemake/7.18
module load slurm
export PYTHONNOUSERSITE=1

# Run Snakemake using SLURM as cluster backend
snakemake \
    --jobs 100 \
    --cores 100 \
    --resources mem_mb=500000 \
    --latency-wait 60 \
    --rerun-incomplete \
    --scheduler ilp \
    --cluster "sbatch --partition=epyc \
                      --ntasks=1 \
                      --cpus-per-task={threads} \
                      --mem={resources.mem_mb} \
                      --time=172:00:00 \
                      -o logs/cluster/{rule}_{wildcards}.out \
                      -e logs/cluster/{rule}_{wildcards}.err"
