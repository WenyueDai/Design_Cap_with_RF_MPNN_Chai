#!/bin/bash
#SBATCH -J step1_rfdiffusion_matrix
#SBATCH --gres=gpu:1
#SBATCH -p ampere
#SBATCH -A GKAMINSKI-SL2-GPU
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem=16g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wd304@cam.ac.uk
#SBATCH --error=/rds/user/wd304/hpc-work/20241101_rfdiffusion_m7/tasks/%A_%a.err
#SBATCH --output=/rds/user/wd304/hpc-work/20241101_rfdiffusion_m7/tasks/%A_%a.out
#SBATCH --array=1-50  # Corrected array specification

module load anaconda
source /usr/local/software/anaconda/3.2019-10/etc/profile.d/conda.sh
conda activate SE3nv-cuda116

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Running outside of SLURM. Setting TASK_ID manually."
    TASK_ID=1
else
    echo "Using SLURM_ARRAY_TASK_ID."
    TASK_ID=$(($SLURM_ARRAY_TASK_ID))
fi

export HYDRA_FULL_ERROR=1

# Get the task command from the tasks file
task=$(sed -n "${TASK_ID}p" /rds/user/wd304/hpc-work/20241101_rfdiffusion_m7/tasks/tasks.txt)

# Execute the task command
echo "Running task: $task"
eval $task
