import itertools
import pandas as pd
from pathlib import Path
import torch
import shutil
import subprocess
import os

"""In this script use the pyrosetta converted protein structure for submission. The pyrosetta structure 
has been tailored with specially positioned valine and glutamate, and we dont set up a chop, only extension
"""

HPC_USER = "wd304@login-icelake.hpc.cam.ac.uk"
REMOTE_DIR = "/rds/user/wd304/hpc-work/20241031_3ult"
LOCAL_FILE_PATH = "/home/eva/0_bury_charged_pair/2_No_initial_de_novo/2_rfdiffusion_cap/20241022_3ult_inpaint_morehelix333000_15V"
INPUT_PDB_SOURCE = "/home/eva/0_bury_charged_pair/2_No_initial_de_novo/2_rfdiffusion_cap/input/3ult_cleaned_renumbered_monomer_modified2.pdb"
LEEWAY = 15 # How much residue is allow to interact with N/C cap (in contact map), and in inpaint
HELIX_LENGTH = 7
GAP_LENGTH = 3

# ----- Step 1: Setup search space -----
cycle_list = range(1, 6)  # Cycle values from 1 to 5
extend_len_list = range(25, 40)  # Extension length from 25 to 40
# Numbers expressed in terms of PyMol residue numbering (i.e. starts at 1, not 0).
crop_NT_list = range(1, 1)  # N-terminal cropping from 1 (no crop) to 1 (max crop). You crop everything before this residue.
crop_CT_list = range(115, 115)  # C-terminal cropping from 114 (max crop) to 115 (no crop). You crop everything after this residue.

# Create the search space using itertools.product
search_space = list(itertools.product(cycle_list, extend_len_list, crop_CT_list, crop_NT_list))

# Create a DataFrame from the search space
df = pd.DataFrame(search_space, columns=['cycle', 'extend_len', 'crop_CT', 'crop_NT'])

# Update chain_A_len: Adjust based on cropping and extension
df['chain_A_len'] = df.apply(lambda row: row['crop_CT'] - (row['crop_NT'] - 1) + 2 * row['extend_len'], axis=1)

# Generate index_str with zero-padded values
df['index_str'] = df.apply(lambda row: f"cropNT_{str(int(row['crop_NT'])).zfill(3)}_cropCT_{str(int(row['crop_CT'])).zfill(3)}_extendlen_{str(int(row['extend_len'])).zfill(3)}_cycle_{str(int(row['cycle'])).zfill(2)}", axis=1)

# Ensure local folders exist
Path(LOCAL_FILE_PATH).mkdir(parents=True, exist_ok=True)
Path(f"{LOCAL_FILE_PATH}/input").mkdir(parents=True, exist_ok=True)
Path(f"{LOCAL_FILE_PATH}/output").mkdir(parents=True, exist_ok=True)
Path(f"{LOCAL_FILE_PATH}/tasks").mkdir(parents=True, exist_ok=True)

# Add 'scaffold_dir', 'rfd_prefix', and 'input_pdb' columns
df['scaffold_dir'] = df.apply(lambda row: f"{REMOTE_DIR}/output/production_run_2/scaffold_dir/{row['index_str']}", axis=1)
df['rfd_prefix'] = df['index_str'].apply(lambda index: f"{REMOTE_DIR}/output/production_run_2/3ult_{index}")
input_pdb = INPUT_PDB_SOURCE.split("/")[-1]
df['input_pdb'] = f"{REMOTE_DIR}/input/{input_pdb}"  # Input PDB file used for all jobs

# ----- Step 1.1: Copy the PDB file to the generated input folder -----

if not os.path.exists(INPUT_PDB_SOURCE):
    print(f"Error: PDB file not found at {INPUT_PDB_SOURCE}. Please check the path.")
else:
    print(f"Copying PDB file from {INPUT_PDB_SOURCE} to {LOCAL_FILE_PATH}/input/{input_pdb}.")
    shutil.copy(INPUT_PDB_SOURCE, f"{LOCAL_FILE_PATH}/input/{input_pdb}")
    print("PDB file copied successfully.")

# Limit the number of scaffold files
num_scaffolds = int(input("Enter the number of scaffold files to generate: "))
df_limited = df.sample(n=num_scaffolds)

# ----- Step 2: Generate intermediate scaffold files locally and write SBATCH task files -----
def get_scaffoldguided(row):
    chain_len = row['chain_A_len']
    block_adj = torch.zeros((chain_len, chain_len)).float()

    # Define N-terminal and C-terminal blocks
    blocks = {
        "N-terminal": (0, row['extend_len'] + 1 + LEEWAY),
        "C-terminal": (chain_len - (row['extend_len'] + LEEWAY), chain_len)
    }
    
    # Populate adjacency matrix for the defined blocks
    for (start, end) in blocks.values():
        for x, y in itertools.product(range(start, end), repeat=2):
            block_adj[x, y] = 1
    
    for line in block_adj.numpy():
        print(" ".join(map(lambda x: str(int(x)), line)))
    
    # Define the alternating pattern for helix and mask
    pattern = [3] * chain_len
    alternating_pattern = [0, 0, 0, 3, 3, 3]
    
    # Apply alternating pattern to N- and C-terminal blocks
    blocks = {
        "N-terminal": (0, blocks["N-terminal"][1] - GAP_LENGTH),
        "C-terminal": (blocks["C-terminal"][0] + GAP_LENGTH, chain_len)
    }
    for start, end in blocks.values():
        pattern[start:end] = [alternating_pattern[i % len(alternating_pattern)] for i in range(end - start)]
    
    print("Pattern: ", pattern)
    tensor_pattern = torch.tensor(pattern).float()

    # Save tensors for scaffold data
    scaffold_local_dir = row['scaffold_dir'].replace(REMOTE_DIR, LOCAL_FILE_PATH)
    Path(scaffold_local_dir).mkdir(parents=True, exist_ok=True)
    torch.save(tensor_pattern, f'{scaffold_local_dir}/{row["index_str"]}_ss.pt')
    torch.save(block_adj, f'{scaffold_local_dir}/{row["index_str"]}_adj.pt')


    scaffold_local_dir = row['scaffold_dir'].replace(REMOTE_DIR, LOCAL_FILE_PATH)
    Path(scaffold_local_dir).mkdir(parents=True, exist_ok=True)
    torch.save(tensor_pattern, f'{scaffold_local_dir}/{row["index_str"]}_ss.pt')
    torch.save(block_adj, f'{scaffold_local_dir}/{row["index_str"]}_adj.pt')

def write_task(row):
    command = ("/home/wd304/.conda/envs/SE3nv-cuda116/bin/python "
               "/rds/user/wd304/hpc-work/RFdiffusion/scripts/run_inference.py "
              f"inference.output_prefix={row['rfd_prefix']} "
              f"inference.input_pdb={row['input_pdb']} "
              f"contigmap.contigs=[{row['extend_len']}/A{row['crop_NT']}-{row['crop_CT']}/{row['extend_len']}] "
               "inference.write_trajectory=False "
               "inference.num_designs=10 "
               "scaffoldguided.scaffoldguided=True "
               "scaffoldguided.target_pdb=False "
               "scaffoldguided.systematic=True "
              f"scaffoldguided.scaffold_dir={row['scaffold_dir']}")

    with open(f'{LOCAL_FILE_PATH}/tasks/tasks.txt', "a") as file:
        file.write(f"{command}\n")

for _, df_row in df_limited.iterrows():
    get_scaffoldguided(df_row)
    write_task(df_row)

# ----- Step 4: Generate SBATCH script locally -----
# SBATCH header with the correct job array specification
sbatch_header = f"""#!/bin/bash
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
#SBATCH --error={REMOTE_DIR}/tasks/%A_%a.err
#SBATCH --output={REMOTE_DIR}/tasks/%A_%a.out
#SBATCH --array=1-{num_scaffolds}  # Corrected array specification

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
task=$(sed -n "${{TASK_ID}}p" {REMOTE_DIR}/tasks/tasks.txt)

# Execute the task command
echo "Running task: $task"
eval $task
"""

# Write the SBATCH script to a file
sbatch_file = f"{LOCAL_FILE_PATH}/tasks/digs_array_job.sh"
with open(sbatch_file, "w") as file:
    file.write(sbatch_header)

print(f"SBATCH script '{sbatch_file}' has been created.")

# ----- Step 5: Upload to the HPC -----
subprocess.run(f"scp -r {LOCAL_FILE_PATH} {HPC_USER}:{REMOTE_DIR}", shell=True)
print(f"SLURM script copied to {REMOTE_DIR}.")

subprocess.run(f"ssh {HPC_USER} 'cd {REMOTE_DIR}/tasks/ && sbatch digs_array_job.sh'", shell=True)
print("Submitted digs_array_job.sh to HPC.")
