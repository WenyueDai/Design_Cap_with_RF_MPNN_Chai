import os
import re
import shutil
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import pymol2
from glob import glob

# Set global font size for plots
plt.rcParams.update({'font.size': 32})

# Configuration
cap_type = "N cap"  # Choose either "N cap" or "C cap"
LEEWAY = 14
MIN_INTERFACE_BSA = 100  # Minimum Interface BSA to filter caps
MIN_HELIX = 80 # Minimum helix percentage

# Paths
input_folder = "/home/eva/0_bury_charged_pair/5_Pipeline_r1_extenlen_25_40/step3_chai_m8/N_cap/analysis_localrmsd/top_files"
output_folder = f"{input_folder}/bsa"
filter_folder = f"{output_folder}/0_filter"
results_file = f"{output_folder}/results.txt"
bsa_histogram_file = f"{output_folder}/bsa_histogram.png"
helix_percentage_file = f"{output_folder}/helix_percentage.png"

# Ensure output folders exist
os.makedirs(output_folder, exist_ok=True)
os.makedirs(filter_folder, exist_ok=True)

# Functions for various criteria

# 1. Function to convert CIF to PDB
def convert_cif_to_pdb(cif_path):
    pdb_path = cif_path.replace(".cif", ".pdb")
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(cif_path, "structure")
        pymol.cmd.save(pdb_path)
    return pdb_path

# 2. Function to get all PDB files, converting CIF files if necessary
def get_pdb_files(input_folder):
    pdb_files = []
    for file in glob(f"{input_folder}/*"):
        if file.endswith(".cif"):
            pdb_files.append(convert_cif_to_pdb(file))
        elif file.endswith(".pdb"):
            pdb_files.append(file)
    return pdb_files

# 3. Function to extract extend_number from filename
def get_extend_number_from_filename(filename):
    match = re.search(r"extendlen_(\d+)", filename)
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"Extend number not found in filename: {filename}")

# 4. Add Hydrogens using PyMOL
def add_hydrogens(pdb_path):
    new_pdb_path_list = []
    for pdb_file in pdb_path:
        new_pdb_path = pdb_file.replace(".pdb", "_with_h.pdb")
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(pdb_file, "structure")
            pymol.cmd.h_add("all")
            pymol.cmd.save(new_pdb_path)
            new_pdb_path_list.append(new_pdb_path)
    return new_pdb_path_list

# 7. Function to calculate Interface Buried Surface Area
def calculate_interface_bsa(pdb_path):
    filename = os.path.basename(pdb_path)
    extend_number = get_extend_number_from_filename(filename)

    try:
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(pdb_path, "structure")

            # Define cap and corresponding body regions
            if cap_type == "N cap":
                cap_residues = f"1-{extend_number}"
                body_start = extend_number + 1
                body_residues = f"{body_start}-{body_start + extend_number - 1}"
            elif cap_type == "C cap":
                cap_residues = f"{extend_number + 1}-last"
                body_start = extend_number - (extend_number - 1)
                body_residues = f"{body_start - extend_number + 1}-{body_start}"
            else:
                print("Cap type not defined.")
                return 0, 0, 0

            # Verify the selections are non-empty
            pymol.cmd.select("cap_region", f"resi {cap_residues}")
            pymol.cmd.select("body_region", f"resi {body_residues}")
            if pymol.cmd.count_atoms("cap_region") == 0 or pymol.cmd.count_atoms("body_region") == 0:
                print(f"Warning: Empty selection in cap or body region for {pdb_path}")
                return 0, 0, 0

            # Calculate surface areas
            pymol.cmd.create("cap_only", "cap_region")
            pymol.cmd.create("body_only", "body_region")
            pymol.cmd.create("combined", "cap_region or body_region")
            pymol.cmd.get_area("cap_only", load_b=1)
            pymol.cmd.get_area("body_only", load_b=1)
            pymol.cmd.get_area("combined", load_b=1)

            # Calculate buried surface area (BSA)
            cap_area = pymol.cmd.get_area("cap_only")
            body_area = pymol.cmd.get_area("body_only")
            combined_area = pymol.cmd.get_area("combined")
            bsa = (cap_area + body_area) - combined_area

            pymol.cmd.delete("all")
        return cap_area, body_area, bsa
    except Exception as e:
        print(f"Error calculating interface BSA for {pdb_path}: {e}")
        return 0, 0, 0

import os
import pymol2

def calculate_helix(pdb_path):
    # Extract the filename and determine the extend number
    filename = os.path.basename(pdb_path)
    extend_number = get_extend_number_from_filename(filename)

    if extend_number is None or extend_number <= 0:
        print("Invalid or missing extend number.")
        return 0

    try:
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(pdb_path, "structure")

            # Define cap and corresponding body regions
            if cap_type == "N cap":
                cap_residues = f"1-{extend_number}"
            elif cap_type == "C cap":
                cap_residues = f"{extend_number + 1}-last"
            else:
                print("Invalid cap type. Choose 'N cap' or 'C cap'.")
                return 0

            # Verify the selections are non-empty
            pymol.cmd.select("cap_region", f"resi {cap_residues}")

            # Perform DSSP secondary structure assignment
            pymol.cmd.dss()

            # Calculate helix percentage
            total_residues = pymol.cmd.count_atoms(f"resi {cap_residues} and name CA")
            helix_residues = pymol.cmd.count_atoms(f"resi {cap_residues} and ss H and name CA")

            if total_residues == 0:
                raise ValueError("No residues found in the cap region.")

            helix_percentage = (helix_residues / total_residues) * 100

            # Cleanup
            pymol.cmd.delete("all")

        return helix_percentage

    except Exception as e:
        print(f"An error occurred: {e}")
        return 0

# 9. Filter and copy PDBs with updated criteria
def filter_and_copy_pdbs(pdb_files):
    data = []
    for pdb_path in pdb_files:
        cap_bsa, body_bsa, bsa = calculate_interface_bsa(pdb_path)
        helix_percentage = calculate_helix(pdb_path)

        data.append([
            os.path.basename(pdb_path), cap_bsa, body_bsa, bsa, helix_percentage
        ])

        # Filtering criteria based on hydrogen bonds, close atoms, and interface BSA/complementarity
        if bsa >= MIN_INTERFACE_BSA and helix_percentage >= MIN_HELIX:
            shutil.copy(pdb_path, f"{filter_folder}")

    df = pd.DataFrame(data, columns=[
        "PDB_File", "Cap_BSA", "Body_BSA", "bsa", "helix_percentage"
    ])
    df.to_csv(results_file, index=False, sep='\t')
    return df

# 10. Plotting functions for histograms
def plot_histograms(df):
    plt.figure(figsize=(10, 10))
    n, bins, patches = plt.hist(df["bsa"], bins=10, alpha=0.7, color='orange')
    plt.xlabel("Interface BSA")
    plt.ylabel("Frequency")
    for i in range(len(n)):
        plt.text(bins[i] + (bins[i+1] - bins[i]) / 2, n[i], str(int(n[i])), 
                 ha='center', va='bottom', fontsize=25)
    plt.title("Interface BSA between Cap and Body")
    plt.savefig(bsa_histogram_file)
    plt.show()
    
    plt.figure(figsize=(10,10))
    n, bins, patches = plt.hist(df["helix_percentage"], bins=5, alpha=0.7, color='orange')
    plt.xlabel("helix percentage")
    plt.ylabel("Frequency")
    for i in range(len(n)):
        plt.text(bins[i] + (bins[i+1] - bins[i]) / 2, n[i], str(int(n[i])), 
                 ha='center', va='bottom', fontsize=25)
    plt.title("Helix percentage")
    plt.savefig(helix_percentage_file)
    plt.show()

# Place this at the end of your script after all functions have run
def cleanup(input_folder):
    # Remove *with_h* files from the input folder
    for file in glob(os.path.join(input_folder, '*with_h*')):
        try:
            os.remove(file)
            print(f"Deleted file: {file}")
        except Exception as e:
            print(f"Error deleting file {file}: {e}")

# Main script to run everything
if __name__ == "__main__":
    pdb_files = get_pdb_files(input_folder)
    pdb_files = add_hydrogens(pdb_files)
    df = filter_and_copy_pdbs(pdb_files)
    plot_histograms(df)
    cleanup(input_folder)
