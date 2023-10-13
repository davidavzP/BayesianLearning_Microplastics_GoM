#!/bin/sh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition comp06
#SBATCH --time=6:00:00

module purge
module --ignore-cache load python/3.10-anaconda
source /share/apps/bin/conda-3.10.sh
conda activate py3_parcels_mpi

python rechunk.py data/output/ data/output/processed/GOM_rt_2021_12.0-7.0_at_180.0.zarr
