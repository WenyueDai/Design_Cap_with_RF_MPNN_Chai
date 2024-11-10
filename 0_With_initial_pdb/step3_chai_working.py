import os
import subprocess
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
import matplotlib
import re
import shutil
import numpy as np

matplotlib.use('Agg')

# Helper function to run shell commands with retries
def shell_run(command, retries=3):
    """Run a shell command with retry attempts."""
    for attempt in range(retries):
        try:
            subprocess.run(command, shell=True, check=True)
            break
        except subprocess.CalledProcessError as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            if attempt == retries - 1:
                raise e

# Function to clean the FASTA headers by keeping only the sample ID
def clean_fasta_headers(input_fasta, output_fasta):
    """Cleans the FASTA headers by formatting them as 'protein|name=sample-X'."""
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            header_part = record.id.split('_')[0]
            cleaned_id = header_part.replace('=', '-')
            formatted_id = f"protein|name={cleaned_id}"
            record.id = formatted_id
            record.description = ""  # Remove the rest of the description
            SeqIO.write(record, output_handle, "fasta")
    return output_fasta

# Write the predict_structure.py script dynamically
def write_predict_structure_script(local_folder):
    """Dynamically generate the predict_structure.py script for HPC."""
    predict_structure_content = f"""
from pathlib import Path
import numpy as np
import torch
from chai_lab.chai1 import run_inference
import sys
import re
import os

# Get the input FASTA file path from the command line
if len(sys.argv) < 2:
    raise ValueError("No FASTA file provided. Usage: python predict_structure.py <fasta_file>")

fasta_file = sys.argv[1]

fasta_basename = os.path.basename(fasta_file)

match = re.search(r"active_site_(\\d+)_(\\d+).", fasta_basename)
sample_match = re.search(r"_sample(\\d+).", fasta_basename)

num1 = int(match.group(1))
num2 = int(match.group(2))
sample_number = int(sample_match.group(1)) if sample_match else None

base_name = f"active_site_{{num1}}_{{num2}}_sample{{sample_number}}"

# Define the output directory
output_dir = Path(f"./outputs/{{base_name}}")

# Ensure the output directory exists
output_dir.mkdir(parents=True, exist_ok=True)

# Run inference on the input FASTA file
candidates = run_inference(
    fasta_file=Path(fasta_file),
    output_dir=output_dir,
    num_trunk_recycles=3,
    num_diffn_timesteps=200,
    seed=42,
    device=torch.device("cuda:0"),
    use_esm_embeddings=True,
)

# Process results
cif_paths = candidates.cif_paths
scores = [rd.aggregate_score for rd in candidates.ranking_data]

# Load additional scores if needed
scores = np.load(output_dir.joinpath("scores.model_idx_2.npz"))
    """
    
    predict_structure_path = os.path.join(local_folder, "predict_structure.py")
    with open(predict_structure_path, "w") as f:
        f.write(predict_structure_content)
    
    return predict_structure_path

# Write the select_best_prediction.py script dynamically
def write_select_best_prediction_script(local_folder):
    """Dynamically generate the select_best_prediction.py script for finding the best predictions."""
    select_best_prediction_content = f"""
import os
import numpy as np
from pathlib import Path
import shutil
import re
import sys

# Function to find the best prediction based on the aggregate score
def find_best_prediction(fasta_file):
    fasta_basename = os.path.basename(fasta_file)

    match = re.search(r"active_site_(\d+)_(\d+).", fasta_basename)
    sample_match = re.search(r"_sample(\d+).", fasta_basename)

    num1 = int(match.group(1))
    num2 = int(match.group(2))
    sample_number = int(sample_match.group(1)) if sample_match else None

    base_name = f"active_site_{{num1}}_{{num2}}_sample{{sample_number}}"
    output_dir = Path(f"./outputs/{{base_name}}")

    npz_files = list(output_dir.glob("scores.model_idx_*.npz"))
    best_score = None
    best_npz = None

    for npz_file in npz_files:
        data = np.load(npz_file)
        aggregate_score = data['aggregate_score'][0]  # Extract aggregate score
        if best_score is None or aggregate_score > best_score:
            best_score = aggregate_score
            best_npz = npz_file

    if best_npz:
        # Match the corresponding CIF file by replacing the 'scores' part with 'pred'
        cif_index = best_npz.stem.split('_')[-1]  # Extract index from npz file (e.g., 0, 1, 2...)
        best_cif = output_dir / f"pred.model_idx_{{cif_index}}.cif"
        return best_npz, best_cif, best_score
    return None, None, None

# Function to move and rename the best prediction files
def move_best_files(fasta_file, best_npz, best_cif, best_score):
    fasta_basename = os.path.basename(fasta_file)

    # Extract the current folder name (e.g., active_site_x_x_samplex)
    match = re.search(r"active_site_(\d+)_(\d+).", fasta_basename)
    sample_match = re.search(r"_sample(\d+).", fasta_basename)

    num1 = int(match.group(1))
    num2 = int(match.group(2))
    sample_number = int(sample_match.group(1)) if sample_match else None

    # The current folder name (base_name) to rename the files
    base_name = f"active_site_{{num1}}_{{num2}}_sample{{sample_number}}"
    output_dir = Path(f"./outputs/{{base_name}}")

    # Define the destination folder for best predictions
    best_predictions_folder = output_dir.parent / "best_prediction"
    best_predictions_folder.mkdir(exist_ok=True)

    # Rename the files with the folder name before moving them
    new_npz_name = f"{{base_name}}.npz"
    new_cif_name = f"{{base_name}}.cif"

    if best_npz.exists() and best_cif.exists():
        # Move the files to the destination folder
        shutil.move(str(best_npz), best_predictions_folder / new_npz_name)
        shutil.move(str(best_cif), best_predictions_folder / new_cif_name)
        print(f"Moved {{new_npz_name}} and {{new_cif_name}} to {{best_predictions_folder}}")
        
        # Write the summary file with best PDB name and aggregate score
        write_best_prediction_summary(base_name, best_predictions_folder, new_cif_name, best_score)
    else:
        print(f"Error: Could not find {{best_npz}} or {{best_cif}}")

# Function to write the best PDB name and aggregate score to a summary text file
def write_best_prediction_summary(base_name, destination_folder, best_pdb_name, best_score):
    summary_file_path = destination_folder / f"best_prediction_summary.txt"
    
    with open(summary_file_path, "a") as summary_file:
        summary_file.write(f"{{best_pdb_name}}: {{best_score }}")
    
    print(f"Summary written to {{summary_file_path}}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python select_best_prediction.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    best_npz, best_cif, best_score = find_best_prediction(fasta_file)
    
    if best_npz and best_cif and best_score is not None:
        move_best_files(fasta_file, best_npz, best_cif, best_score)
    else:
        print("No valid predictions found in", fasta_file)
    """
    
    select_best_prediction_path = os.path.join(local_folder, "select_best_prediction.py")
    with open(select_best_prediction_path, "w") as f:
        f.write(select_best_prediction_content)
    
    return select_best_prediction_path

# Function to process all FASTA files and submit jobs to HPC
def process_fasta_files(fasta_files, job_name):
    local_folder = os.path.join(os.getcwd(), job_name)
    hpc_user = "wd304@login-icelake.hpc.cam.ac.uk"
    hpc_remote_dir = f"/rds/user/wd304/hpc-work/{job_name}"

    os.makedirs(local_folder, exist_ok=True)

    cleaned_fasta_files = []
    hpc_fasta_files = []
    for fasta_file_path in fasta_files:
        clean_fasta_path = os.path.join(local_folder, f"cleaned_{os.path.basename(fasta_file_path)}")
        clean_fasta_headers(fasta_file_path, clean_fasta_path)  # Clean the headers
        cleaned_fasta_files.append(clean_fasta_path)
        
        # Write predict_structure.py into the local folder
        write_predict_structure_script(local_folder)
        
        # Write select_best_prediction.py into the local folder
        write_select_best_prediction_script(local_folder)

        # HPC path for FASTA
        hpc_fasta_path = f"{hpc_remote_dir}/cleaned_{os.path.basename(fasta_file_path)}"
        hpc_fasta_files.append(hpc_fasta_path)

    # Generate SLURM script with HPC paths
    slurm_script_content, slurm_script_path = create_slurm_script(hpc_remote_dir, job_name, hpc_fasta_files)
    
    with open(os.path.join(local_folder, f"{job_name}.sh"), "w") as f:
        f.write(slurm_script_content)

    # Upload the folder to HPC and submit job
    upload_to_hpc(local_folder, hpc_user, hpc_remote_dir)
    submit_job_array_to_hpc(hpc_user, hpc_remote_dir, job_name)

# SLURM and upload functions
def create_slurm_script(hpc_remote_dir, job_name, hpc_fasta_files):
    fasta_files_str = ' '.join([f"'{f}'" for f in hpc_fasta_files])  # Add HPC path to the files

    slurm_script_content = f"""#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH -p ampere
#SBATCH -A GKAMINSKI-SL2-GPU
#SBATCH -t 20:00:00
#SBATCH --nodes=1
#SBATCH --array=0-{len(hpc_fasta_files) - 1}%10  # Process 2 FASTA files at once
#SBATCH -o {hpc_remote_dir}/job_output_%A_%a.log
#SBATCH -e {hpc_remote_dir}/job_error_%A_%a.log

module add anaconda
source /usr/local/software/anaconda/3.2019-10/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate chai_py310

fasta_files=({fasta_files_str})
fasta_file=${{fasta_files[$SLURM_ARRAY_TASK_ID]}}

/home/wd304/.conda/envs/chai_py310/bin/python {hpc_remote_dir}/predict_structure.py "$fasta_file"
/home/wd304/.conda/envs/chai_py310/bin/python {hpc_remote_dir}/select_best_prediction.py "$fasta_file"
"""
    return slurm_script_content, os.path.join(hpc_remote_dir, f"{job_name}.sh")

def submit_job_array_to_hpc(hpc_user, hpc_remote_dir, job_name):
    shell_run(f"ssh {hpc_user} 'cd {hpc_remote_dir} && sbatch {job_name}.sh'")

def upload_to_hpc(local_folder, hpc_user, hpc_remote_dir):
    response = input("Do you want to upload the generated files to the HPC? (y/n): ").lower()
    if response == 'y':
        shell_run(f"scp -r {local_folder} {hpc_user}:{hpc_remote_dir}")
        print(f"Files uploaded to {hpc_user}:{hpc_remote_dir}")
    else:
        print("Upload cancelled.")

# Main function to handle job submissions
def create_submission_scripts(fasta_dir):
    """Create and submit SLURM scripts for a list of FASTA files."""
    folder_name = os.path.basename(os.path.normpath(fasta_dir))
    job_name = f"chai_{folder_name}_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}"

    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith(".fasta")]

    if not fasta_files:
        print("No FASTA files found in the directory.")
        return

    # Process the FASTA files using the job array
    process_fasta_files(fasta_files, job_name)

if __name__ == "__main__":
    # Input the folder path containing FASTA files
    fasta_dir = input("Enter the local folder path containing FASTA files: ")

    # Check if the folder exists
    if not os.path.isdir(fasta_dir):
        print(f"Error: Folder '{fasta_dir}' does not exist.")
        exit(1)

    # Start the job submission process
    create_submission_scripts(fasta_dir)
