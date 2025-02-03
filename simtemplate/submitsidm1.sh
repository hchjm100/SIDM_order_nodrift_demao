#!/bin/bash -l

#SBATCH --ntasks=64
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --time=7-00:00:00     # 7 day and 15 minutes
#SBATCH --mail-user=dkong012@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="sidm10id18dmo"
#SBATCH -p batch # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

# Print current date
date

# Load samtools
module load gadget2
module load hdf5/1.12.1
module load fftw/2.1.5
#module load mpich/4.0.1
module load openmpi/4.1.2_slurm-21.08.5

# Concatenate BAMs
mpirun -np 64 /rhome/dkong012/bigdata/hisim/dmo_run_id18_sidm_10/codeSIDM1/GadgetSIDM1 sidm1.param
#samtools cat -h header.sam -o out.bam in1.bam in2.bam

# Print name of node
hostname
