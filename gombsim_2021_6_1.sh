#!/bin/sh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition comp72
#SBATCH --time=20:00:00

module purge
module --ignore-cache load python/3.10-anaconda
source /share/apps/bin/conda-3.10.sh
conda activate py3_parcels_mpi

python gom_b_sim_v2.py data/output_v2/raw4/ 2021 6 1 180

python rechunk.py data/output_v2/raw4/ data/output_v2/rechunked/GOM_rt_2021_6.0-1.0_at_180.0.zarr


