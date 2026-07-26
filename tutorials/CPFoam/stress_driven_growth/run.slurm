#!/bin/bash
#SBATCH --job-name=cpfoam_growth
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=04:00:00
#SBATCH --output=log.slurm-%j.out
#SBATCH --error=log.slurm-%j.err

module purge
module load apps/gcc/openfoam/v2506
source "$foamDotFile"

cd "$SLURM_SUBMIT_DIR" || exit 1

# Sync decomposeParDict with actual task count
sed -i "s/^numberOfSubdomains.*/numberOfSubdomains ${SLURM_NTASKS};/" system/decomposeParDict

./Allrun parallel
