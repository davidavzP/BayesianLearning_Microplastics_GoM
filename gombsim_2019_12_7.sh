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

python gom_b_sim_v2.py data/output_v2/raw1/ 2019 12 7 180

python rechunk.py data/output_v2/raw1/ data/output_v2/rechunked/GOM_rt_2019_12.0-7.0_at_180.0.zarr

# hostfile=/scratch/${SLURM_JOB_ID}/machinefile_${SLURM_JOB_ID}

# module purge
# module load gcc/11.2.1 mkl/19.0.5 python/3.10-anaconda
# source /share/apps/bin/conda-3.10.sh
# conda activate py3_parcels_mpi

# module load mvapich2/2.3.7 cuda/11.7
# mpiexec -ppn 16 -hostfile $hostfile python gom_b_sim_no_rerun.py data/output/raw1/ 2019 12 7 180

# python rechunk.py data/output/raw1/GOM_rt_2019_12.0-7.0_at_180.0/ data/output/processed/GOM_rt_2019_12.0-7.0_at_180.0.zarr




