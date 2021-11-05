#!/bin/bash -login
#SBATCH -p high                # partition, or queue, to assign to
#SBATCH -J gse_controls         # name for job
#SBATCH -N 1                   # one "node", or computer
#SBATCH -n 1                   # one task for this node
#SBATCH -c 1                   # eight cores per task
#SBATCH -t 30:00:00            # ask for no more than 30 minutes
#SBATCH --mem=16Gb             # ask for no more than 10 GB of memory

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate your desired conda environment
conda activate singlecell-dev

# go to the directory you ran 'sbatch' in, OR just hardcode it...

# fail on weird errors
set -o nounset
set -o errexit
set -x

# Run r scripts
#Rscript $HOME/ets04-2/R/grouped-sccatch.R # rna counts
Rscript $HOME/2021-geo_datasets-gse149758/R/catch.r # rna velocity

# print out various information about the job
env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch