#!/bin/sh
#SBATCH --mem=30G
#SBATCH --time=01-00:00:00

srun R CMD BATCH hierarchical_model.R
