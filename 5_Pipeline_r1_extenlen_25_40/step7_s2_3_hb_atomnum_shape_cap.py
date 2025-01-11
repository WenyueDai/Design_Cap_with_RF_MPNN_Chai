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
MIN_HBONDS = 8  # Minimum number of hydrogen bonds to keep the file
MIN_CLOSE_ATOMS = 180
MIN_INTERFACE_BSA = 140  # Minimum Interface BSA to filter caps
MIN_HYDROPHOBIC_CONTACTS = 20  # Threshold for hydrophobic contacts

DISTANCE_CUTOFF = 3.5  # Distance cutoff for hydrogen bonds in Å
ANGLE_CUTOFF = 135.0  # Angle cutoff for hydrogen bonds in degrees
CLOSE_ATOM_DISTANCE = 3

# Paths
input_folder = "/home/eva/0_bury_charged_pair/5_Pipeline_r1_extenlen_25_40/step3_chai_m8/N_cap/analysis_localrmsd/top_files"
output_folder = f"{input_folder}/hbond"
cpptraj_output_folder = f"{output_folder}/cpptraj_outputs"
filter_folder = f"{output_folder}/0_filter"
results_file = f"{output_folder}/results.txt"
scatter_plot_file = f"{output_folder}/scatter_plot.png"
hbond_histogram_file = f"{output_folder}/hbond_histogram.png"
close_atom_histogram_file = f"{output_folder}/close_atom_histogram.png"
bsa_histogram_file = f"{output_folder}/bsa_histogram.png"
complementarity_histogram_file = f"{output_folder}/complementarity_histogram.png"

# Ensure output folders exist
os.makedirs(cpptraj_output_folder, exist_ok=True)
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

# 5. Function to calculate hydrogen bonds with cpptraj, including backbone-only calculation
import os
import re
import subprocess

# 5. Function to calculate hydrogen bonds with cpptraj, including backbone-only calculation
def calculate_hbonds_with_cpptraj(pdb_path):
    filename = os.path.basename(pdb_path)
    extend_number = get_extend_number_from_filename(filename)

    # Define cap and leeway regions based on cap type
    if cap_type == "N cap":
        cap_residues = f"1-{extend_number}"
        leeway_residues = f"{extend_number + 1}-{extend_number + LEEWAY}"
    elif cap_type == "C cap":
        cap_residues = f"{extend_number + 1}-last"
        leeway_residues = f"{extend_number - LEEWAY + 1}-{extend_number}"
    else:
        print("Cap type not defined. Skipping hydrogen bond analysis.")
        return 0, 0, 0

    def run_cpptraj(input_commands, output_file):
        cpptraj_input_file = "cpptraj_temp.in"
        with open(cpptraj_input_file, "w", newline="\n") as file:
            file.write("\n".join(input_commands))
        try:
            subprocess.run(["cpptraj", "-i", cpptraj_input_file], capture_output=True, text=True, check=True)
            os.remove(cpptraj_input_file)  # Clean up temporary input file
            
            # Initialize hbond count by summing values from each frame in the output
            hbond_count = 0
            with open(output_file, "r") as output:
                for line in output:
                    # Look for lines with data (e.g., "1           14")
                    match = re.search(r'^\s*\d+\s+(\d+)', line)
                    if match:
                        hbond_count += int(match.group(1))  # Add the hydrogen bond count from this frame

            return hbond_count
        except Exception as e:
            print(f"Error running cpptraj for {pdb_path}: {e}")
            if os.path.exists(cpptraj_input_file):
                os.remove(cpptraj_input_file)  # Ensure temp file is deleted on error
            return 0

    # Define output files for each type of hydrogen bond calculation
    combined_output = f"{cpptraj_output_folder}/{filename}_combined_hbonds.dat"
    backbone_output = f"{cpptraj_output_folder}/{filename}_backbone_hbonds.dat"

    # Define cpptraj commands
    cpptraj_input_combined = [
        f"parm {pdb_path}",
        f"trajin {pdb_path}",
        f"hbond HB out {combined_output} :{cap_residues},{leeway_residues} distance {DISTANCE_CUTOFF} angle {ANGLE_CUTOFF}",
        "run", "exit"
    ]
    cpptraj_input_backbone = [
        f"parm {pdb_path}",
        f"trajin {pdb_path}",
        f"hbond Backbone :{cap_residues},{leeway_residues}@C,O,N,H out {backbone_output} distance {DISTANCE_CUTOFF} angle {ANGLE_CUTOFF}",
        "run", "exit"
    ]

    # Run cpptraj commands and get bond counts
    total_hbonds = run_cpptraj(cpptraj_input_combined, combined_output)
    backbone_hbonds = run_cpptraj(cpptraj_input_backbone, backbone_output)

    # Calculate non-backbone hydrogen bonds
    non_backbone_hbonds = max(0, total_hbonds - backbone_hbonds)

    # Print results for verification
    print(f"{filename}: Total H-Bonds = {total_hbonds}, Backbone H-Bonds = {backbone_hbonds}, Non-Backbone H-Bonds = {non_backbone_hbonds}")

    return non_backbone_hbonds, total_hbonds, backbone_hbonds

# 7. Function to calculate Interface Buried Surface Area
import os
import pymol2

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


# 8. Function to calculate close atoms between cap and leeway using only side-chain atoms in PyMOL
def calculate_close_atoms_pymol(pdb_path):
    filename = os.path.basename(pdb_path)
    extend_number = get_extend_number_from_filename(filename)

    # Define cap and leeway regions based on cap type
    if cap_type == "N cap":
        cap_residues = f"1-{extend_number}"
        leeway_residues = f"{extend_number + 1}-{extend_number + LEEWAY}"
    elif cap_type == "C cap":
        cap_residues = f"{extend_number + 1}-last"
        leeway_residues = f"{extend_number - LEEWAY + 1}-{extend_number}"
    else:
        print("Cap type not defined.")
        return 0

    try:
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(pdb_path, "structure")
            
            # Select only side-chain atoms in the cap and leeway regions
            pymol.cmd.select("cap_region", f"resi {cap_residues}")
            pymol.cmd.select("leeway_region", f"resi {leeway_residues}")
            
            # Select atoms in cap_region that are within the cutoff distance of leeway_region
            pymol.cmd.select("close_atoms", f"br. cap_region within {CLOSE_ATOM_DISTANCE} of leeway_region")
            close_atom_count = pymol.cmd.count_atoms("close_atoms")
            pymol.cmd.delete("all")
            
        return close_atom_count
    except Exception as e:
        print(f"Error calculating close atoms for {pdb_path}: {e}")
        return 0

# 9. Filter and copy PDBs with updated criteria
def filter_and_copy_pdbs(pdb_files):
    data = []
    for pdb_path in pdb_files:
        non_mainchain_hbonds, total_hbonds, mainchain_hbonds = calculate_hbonds_with_cpptraj(pdb_path)
        num_close_atoms = calculate_close_atoms_pymol(pdb_path)
        cap_bsa, body_bsa, bsa = calculate_interface_bsa(pdb_path)

        data.append([
            os.path.basename(pdb_path), non_mainchain_hbonds, total_hbonds, mainchain_hbonds,
            num_close_atoms, cap_bsa, body_bsa, bsa
        ])

        # Filtering criteria based on hydrogen bonds, close atoms, and interface BSA/complementarity
        if non_mainchain_hbonds >= MIN_HBONDS and num_close_atoms >= MIN_CLOSE_ATOMS and \
           bsa >= MIN_INTERFACE_BSA:
            shutil.copy(pdb_path, f"{filter_folder}")

    df = pd.DataFrame(data, columns=[
        "PDB_File", "Non_MainChain_H-Bonds", "Total_H-Bonds", "MainChain_H-Bonds",
        "Close_Atoms", "Cap_BSA", "Body_BSA", "bsa"
    ])
    df.to_csv(results_file, index=False, sep='\t')
    return df

# 10. Plotting functions for histograms
def plot_histograms(df):
    plt.figure(figsize=(10, 10))
    n, bins, patches = plt.hist(df["Non_MainChain_H-Bonds"], bins=10, alpha=0.7, color='orange')
    plt.xlabel("Non Main-Chain H-Bonds")
    plt.ylabel("Frequency")
    for i in range(len(n)):
        plt.text(bins[i] + (bins[i+1] - bins[i]) / 2, n[i], str(int(n[i])), 
                 ha='center', va='bottom', fontsize=25)
    plt.title("Non Main-Chain H-Bonds between Cap and Leeway")
    plt.savefig(hbond_histogram_file)
    plt.show()

    plt.figure(figsize=(10, 10))
    n, bins, patches = plt.hist(df["Close_Atoms"], bins=10, alpha=0.7, color='orange')
    plt.xlabel("Close Atom Count")
    plt.ylabel("Frequency")
    for i in range(len(n)):
        plt.text(bins[i] + (bins[i+1] - bins[i]) / 2, n[i], str(int(n[i])), 
                 ha='center', va='bottom', fontsize=25)
    plt.title("Close Atoms within 3Å between Cap and Leeway")
    plt.savefig(close_atom_histogram_file)
    plt.show()

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


# Place this at the end of your script after all functions have run
def cleanup(input_folder, cpptraj_output_folder):
    # Remove *with_h* files from the input folder
    for file in glob(os.path.join(input_folder, '*with_h*')):
        try:
            os.remove(file)
            print(f"Deleted file: {file}")
        except Exception as e:
            print(f"Error deleting file {file}: {e}")

    # Remove the entire cpptraj_output_folder
    try:
        shutil.rmtree(cpptraj_output_folder)
        print(f"Deleted folder and its contents: {cpptraj_output_folder}")
    except Exception as e:
        print(f"Error deleting folder {cpptraj_output_folder}: {e}")

# Main script to run everything
if __name__ == "__main__":
    pdb_files = get_pdb_files(input_folder)
    pdb_files = add_hydrogens(pdb_files)
    df = filter_and_copy_pdbs(pdb_files)
    plot_histograms(df)
    cleanup(input_folder, cpptraj_output_folder)
