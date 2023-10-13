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

python gom_b_sim_v2.py data/output_v2/raw3/ 2020 12 7 180

python rechunk.py data/output_v2/raw3/ data/output_v2/rechunked/GOM_rt_2020_12.0-7.0_at_180.0.zarr

# python gom_b_sim.py data/output/raw3/ 2020 12 7 180

# python rechunk.py data/output/raw3/ data/output/processed/GOM_rt_2020_12.0-7.0_at_180.0.zarr


