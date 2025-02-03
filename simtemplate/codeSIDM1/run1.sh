#!/bin/bash -l
#
#SBATCH --output=GadgetSIDM1_run1.out     # Standard output and error log
#
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=64              # Total MPI tasks (matching mpirun -np #)
#SBATCH --cpus-per-task=1       # Number of CPU cores per MPI task
#SBATCH --mem=30G                # Total memory per node (adjust as needed)
#SBATCH --time=10-00:00:00       # Time limit (D-HH:MM:SS) this case: 10 days
#SBATCH --mail-user=hsang012@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="1order1_demao"
#SBATCH -p batch                 # Partition (change to epyc, intel, batch, etc.)

module load fftw/2.1.5

# 2) Print some diagnostic info
echo "Starting GadgetSIDM1 job on ${HOSTNAME}"
echo "Date: $(date)"
echo "Directory: $(pwd)"

# 3) Run your program
# If the executable and .param file are in the current directory:
mpirun -np 64 ./GadgetSIDM1 run1.param

# 4) Print completion time
echo "Finished at: $(date)"

