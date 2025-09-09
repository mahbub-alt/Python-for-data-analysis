#!/bin/bash
#SBATCH --job-name=serial_job_test    # Job name
#SBATCH --partition=water           # Partition Name (Required)
#SBATCH --mail-type=END       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mahbubenv@ku.edu      # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=64g 
#SBATCH --time=10-6:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log


# module load python
module load conda
echo "Running 30p_python script"
# conda activate /kuhpc/work/water/m089r172/conda/envs/new_envv 
conda activate aq_env

python /kuhpc/home/m089r172/my_first_exmpl/to_run_NSGA-II.py 