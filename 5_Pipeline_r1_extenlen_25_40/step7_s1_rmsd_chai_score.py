import os
import shutil
import re
import subprocess
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.exceptions import SelectionError
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set global font size for plots
plt.rcParams.update({'font.size': 32})

cap_type = "N cap"  # Choose either "N cap" or "C cap"

# Define the paths
chai_folder = "/home/eva/0_bury_charged_pair/5_Pipeline_r1_extenlen_25_40/step3_chai_m8/N_cap"
rdfiffusion_folder = "/home/eva/0_bury_charged_pair/5_Pipeline_r1_extenlen_25_40/step1_rfdiffusion_m8/1_cap_calculation/N_cap_add_3ultback"
analysis_folder = os.path.join(chai_folder, "analysis_localrmsd")
rmsd_output_file = os.path.join(analysis_folder, "rmsd_results.txt")
scatter_plot_file = os.path.join(analysis_folder, "rmsd_vs_chai_score.png")
rmsd_vs_chai_output_file = os.path.join(analysis_folder, "rmsd_vs_chai_score.txt")
rmsd_vs_extendlen_file = os.path.join(analysis_folder, "rmsd_vs_extendlen.png")
summary_file = os.path.join(chai_folder, "best_prediction_summary.txt")

# Ensure that the analysis folder exists
os.makedirs(analysis_folder, exist_ok=True)

rmsd_data = []
chai_scores = []
extend_lengths = []
sorted_data = []

# Function to calculate global RMSD for Cα atoms only
def calculate_global_rmsd(structure1, structure2, extend_number):
    try:
        u1 = mda.Universe(structure1)
        u2 = mda.Universe(structure2)
        
        # Get the number of residues in the structure
        num_residues_u1 = len(u1.select_atoms("name CA").residues)
        num_residues_u2 = len(u2.select_atoms("name CA").residues)
        
        LEEWAY = 14

        # Select residues based on cap type
        if cap_type == "N cap":
            u1_selection = u1.select_atoms(f"name CA and resid 1:{extend_number + LEEWAY}")
            u2_selection = u2.select_atoms(f"name CA and resid 1:{extend_number + LEEWAY}")
        elif cap_type == "C cap":
            start_res_u1 = max(num_residues_u1 - (extend_number + LEEWAY), 1)
            start_res_u2 = max(num_residues_u2 - (extend_number + LEEWAY), 1)
            u1_selection = u1.select_atoms(f"name CA and resid {start_res_u1}:{num_residues_u1}")
            u2_selection = u2.select_atoms(f"name CA and resid {start_res_u2}:{num_residues_u2}")
        else:  # Full structure comparison if cap_type is None or not specified
            u1_selection = u1.select_atoms("name CA")
            u2_selection = u2.select_atoms("name CA")
            print("Comparing full structure.")

        # Print the residues being compared
        print(f"Comparing residues for RMSD: {u1_selection.residues} with {u2_selection.residues}")

        alignment = rms.RMSD(u1_selection, u2_selection)
        alignment.run()
        return alignment.results.rmsd[-1, 2]
    except SelectionError as e:
        print(f"SelectionError: {e}")
        return None

# Function to convert .cif to .pdb using PyMOL
def convert_cif_to_pdb(cif_file, output_folder):
    pdb_file = os.path.join(output_folder, os.path.splitext(os.path.basename(cif_file))[0] + ".pdb")
    try:
        pymol_cmd = f"pymol -c -d 'load {cif_file}; save {pdb_file}; quit'"
        subprocess.run(pymol_cmd, shell=True, check=True)
        return pdb_file
    except subprocess.CalledProcessError as e:
        print(f"Error converting {cif_file} to PDB using PyMOL: {e}")
        return None

# Function to extract details from filename
def extract_details_from_filename(filename):
    match = re.search(r"cycle_(\d+)_(\d+)", filename)
    sample_match = re.search(r"_sample(\d+)", filename)
    extend_match = re.search(r"extendlen_(\d+)", filename)
    
    num1 = int(match.group(1))
    num2 = int(match.group(2))
    sample_number = int(sample_match.group(1))
    extend_number = int(extend_match.group(1))

    if cap_type == "N cap":
        active_site_part = f"extendlen_{extend_number}_cycle_{num1}_{num2}_N"
    elif cap_type == "C cap":
        active_site_part = f"extendlen_{extend_number}_cycle_{num1}_{num2}_C"
    else:
        active_site_part = f"extendlen_{extend_number}_cycle_{num1}_{num2}"

    return active_site_part, num1, num2, sample_number, extend_number

# Function to read Chai scores from the summary file
def read_chai_scores(summary_file):
    chai_scores = {}
    with open(summary_file, "r") as f:
        lines = re.findall(r"(extendlen_\d+_cycle_\d+_\d+_sample\d+\.cif):\s+([\d\.]+)", f.read())
        for file_name, score in lines:
            chai_scores[file_name] = float(score)
    return chai_scores

# New function to plot RMSD vs extendlen
def plot_rmsd_vs_extendlen(rmsd_data, extend_lengths):
    plt.figure(figsize=(10, 6))
    plt.scatter(extend_lengths, rmsd_data, alpha=0.7, color='green')
    plt.xlabel("Extendlen")
    plt.ylabel("Global RMSD (Å)")
    plt.title("RMSD vs. Extendlen")
    plt.tight_layout()
    plt.savefig(rmsd_vs_extendlen_file, format="png", dpi=300)
    plt.show()
    print(f"Scatter plot of RMSD vs. Extendlen saved as {rmsd_vs_extendlen_file}.")

# Modified `copy_top_files` function to handle `top_n="all"`
def copy_top_files(sorted_file, source_folder, destination_folder, top_n=30):
    try:
        df_sorted = pd.read_csv(sorted_file, sep='\t')
        df_sorted['global_rmsd'] = pd.to_numeric(df_sorted['global_rmsd'], errors='coerce')
        
        # Drop rows with NaN in 'global_rmsd'
        df_sorted = df_sorted.dropna(subset=['global_rmsd'])
        
        # Determine files to copy
        if top_n == "all":
            top_files = df_sorted['file_name'].tolist()
        else:
            top_files = df_sorted.nsmallest(top_n, 'global_rmsd')['file_name'].tolist()
        
        os.makedirs(destination_folder, exist_ok=True)
        
        # Copy selected files
        for file_name in top_files:
            source_file_path = os.path.join(source_folder, file_name)
            if os.path.exists(source_file_path):
                shutil.copy(source_file_path, destination_folder)
                print(f"Copied {file_name} to {destination_folder}.")
            else:
                print(f"File {file_name} does not exist in {source_folder}.")
    except Exception as e:
        print(f"Error copying top files: {e}")

# Main processing loop (Calculate RMSD, Read Scores, Generate Plots)
chai_scores_dict = read_chai_scores(summary_file)

with open(rmsd_output_file, "w") as rmsd_file, open(rmsd_vs_chai_output_file, "w") as rmsd_vs_chai_file:
    rmsd_file.write("File Name\tGlobal RMSD (Å)\tChai Score\n")
    rmsd_vs_chai_file.write("Extendlen Num\tCycle Num Num\tSample Num\tGlobal RMSD (Å)\tChai Score\n")

    for chai_file in os.listdir(chai_folder):
        if chai_file.endswith(".cif"):
            chai_file_path = os.path.join(chai_folder, chai_file)
            active_site_part, num1, num2, sample_number, extend_number = extract_details_from_filename(chai_file)

            if active_site_part:
                matching_pdb_files = [
                    f for f in os.listdir(rdfiffusion_folder)
                    if active_site_part is not None and re.search(rf"extendlen_0*{extend_number}_cycle_0*{num1}_{num2}", f) and f.endswith(".pdb")
                ]

                if matching_pdb_files:
                    pdb_file_path = os.path.join(rdfiffusion_folder, matching_pdb_files[0])
                    print(f"Comparing: {chai_file} with {pdb_file_path}")

                    # Convert .cif to .pdb for chai_file
                    converted_pdb_path = convert_cif_to_pdb(chai_file_path, chai_folder)
                    if converted_pdb_path:
                        global_rmsd = calculate_global_rmsd(pdb_file_path, converted_pdb_path, extend_number)
                        chai_score = chai_scores_dict.get(chai_file, None)

                        if global_rmsd is not None and chai_score is not None:
                            rmsd_data.append(global_rmsd)
                            chai_scores.append(chai_score)
                            extend_lengths.append(extend_number)

                            rmsd_file.write(f"{chai_file}\t{global_rmsd:.4f}\t{chai_score:.4f}\n")
                            rmsd_vs_chai_file.write(f"{active_site_part}\t{sample_number}\t{global_rmsd:.4f}\t{chai_score:.4f}\n")
                            sorted_data.append((active_site_part, sample_number, chai_file, global_rmsd, chai_score, extend_number))
                else:
                    print(f"No matching file found for {chai_file} in RFdiffusion folder.")
            else:
                print(f"Could not extract details from filename: {chai_file}")

# Plot RMSD vs. Chai score
if rmsd_data and chai_scores:
    plt.figure(figsize=(10, 6))
    plt.scatter(rmsd_data, chai_scores, alpha=0.7, color='b')
    plt.xlabel("Global RMSD (Å)")
    plt.ylabel("Chai Score")
    plt.title("RMSD vs. Chai Score")
    plt.ylim(0.12, 0.2)
    plt.xlim(0, 20)
    plt.tight_layout()
    plt.savefig(scatter_plot_file, format="png", dpi=300)
    plt.show()
    print(f"Scatter plot saved as {scatter_plot_file}.")

# Plot RMSD vs Extendlen
plot_rmsd_vs_extendlen(rmsd_data, extend_lengths)

# Sort and save RMSD and Chai score data
if sorted_data:
    df_sorted = pd.DataFrame(sorted_data, columns=["extend_cycle_num_num", "sample_num", "file_name", "global_rmsd", "chai_score", "extend_length"])
    df_sorted = df_sorted[(df_sorted['global_rmsd'] <= 6) & (df_sorted['chai_score'] >= 0.16)]

    sorted_output_file = os.path.join(analysis_folder, "sorted_rmsd_chai_scores.txt")
    df_sorted.to_csv(sorted_output_file, index=False, sep='\t')
    print(f"Sorted RMSD and Chai score data saved to {sorted_output_file}.")

# Define the folder to copy the top files to
top_files_folder = os.path.join(analysis_folder, "top_files")

# Copy top files (set top_n="all" to copy all filtered files, top_n=10 to copy top 10)
copy_top_files(sorted_output_file, chai_folder, top_files_folder, top_n="all")
