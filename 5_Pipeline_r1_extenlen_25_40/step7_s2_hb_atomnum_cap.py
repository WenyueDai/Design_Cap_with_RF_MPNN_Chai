import os
import re
import shutil
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import pymol2
from glob import glob

# Conda activate AmberTools23

# Set global font size for plots
plt.rcParams.update({'font.size': 32})

# Configuration
cap_type = "N cap"  # Choose either "N cap" or "C cap"
LEEWAY = 14
MIN_HBONDS = 4  # Minimum number of hydrogen bonds to keep the file
MIN_CLOSE_ATOMS = 180

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

# Ensure output folders exist
os.makedirs(cpptraj_output_folder, exist_ok=True)
os.makedirs(output_folder, exist_ok=True)
os.makedirs(filter_folder, exist_ok=True)

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

# 5. Function to calculate hydrogen bonds with cpptraj
def calculate_hbonds_with_cpptraj(pdb_path):
    filename = os.path.basename(pdb_path)
    extend_number = get_extend_number_from_filename(filename)

    if cap_type == "N cap":
        cap_residues = f"1-{extend_number}"
        leeway_residues = f"{extend_number + 1}-{extend_number + LEEWAY}"
    elif cap_type == "C cap":
        cap_residues = f"{extend_number + 1}-last"
        leeway_residues = f"{extend_number - LEEWAY + 1}-{extend_number}"
    else:
        print("Cap type not defined. Skipping hydrogen bond analysis.")
        return 0, 0, 0, 0

    def run_cpptraj(input_commands):
        cpptraj_input_file = "cpptraj_temp.in"
        with open(cpptraj_input_file, "w", newline="\n") as file:
            file.write("\n".join(input_commands))
        try:
            result = subprocess.run(["cpptraj", "-i", cpptraj_input_file], capture_output=True, text=True, check=True)
            os.remove(cpptraj_input_file)  # Clean up temporary input file
            match = re.search(r'(\d+) solute-solute hydrogen bonds.', result.stdout)
            if match:
                return int(match.group(1))  # Extract HB count from output
            else:
                return 0  # Default to 0 if no match is found
        except Exception as e:
            print(f"Error running cpptraj for {pdb_path}: {e}")
            os.remove(cpptraj_input_file)  # Ensure temp file is deleted on error
            return 0

    # Define cpptraj commands
    cpptraj_input_combined = [
        f"parm {pdb_path}", f"trajin {pdb_path}",
        f"hbond HB out {cpptraj_output_folder}/dummy.dat :{cap_residues},{leeway_residues} distance {DISTANCE_CUTOFF} angle {ANGLE_CUTOFF}",
        "run", "exit"]
    cpptraj_input_cap = [
        f"parm {pdb_path}", f"trajin {pdb_path}",
        f"hbond HB out {cpptraj_output_folder}/dummy.dat :{cap_residues} distance {DISTANCE_CUTOFF} angle {ANGLE_CUTOFF}",
        "run", "exit"]
    cpptraj_input_leeway = [
        f"parm {pdb_path}", f"trajin {pdb_path}",
        f"hbond HB out {cpptraj_output_folder}/dummy.dat :{leeway_residues} distance {DISTANCE_CUTOFF} angle {ANGLE_CUTOFF}",
        "run", "exit"]

    # Run cpptraj and capture hbond counts
    hbond_combined = run_cpptraj(cpptraj_input_combined)
    hbond_cap = run_cpptraj(cpptraj_input_cap)
    hbond_leeway = run_cpptraj(cpptraj_input_leeway)

    # Calculate the number of hydrogen bonds between cap and leeway regions
    hbond_between = max(0, hbond_combined - (hbond_cap + hbond_leeway))
    
    return hbond_between, hbond_combined, hbond_cap, hbond_leeway

# 6. Function to calculate close atoms between cap and leeway using PyMOL
def calculate_close_atoms_pymol(pdb_path):
    filename = os.path.basename(pdb_path)
    extend_number = get_extend_number_from_filename(filename)

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
            pymol.cmd.select("cap_region", f"resi {cap_residues}")
            pymol.cmd.select("leeway_region", f"resi {leeway_residues}")
            pymol.cmd.select("close_atoms", f"br. cap_region within {CLOSE_ATOM_DISTANCE} of leeway_region")
            close_atom_count = pymol.cmd.count_atoms("close_atoms")
            pymol.cmd.delete("all")
        return close_atom_count
    except Exception as e:
        print(f"Error calculating close atoms for {pdb_path}: {e}")
        return 0

# 7. Function to process PDB files and save results to DataFrame
def filter_and_copy_pdbs(pdb_files):
    data = []
    for pdb_path in pdb_files:
        hbond_between, hbond_combined, hbond_cap, hbond_leeway = calculate_hbonds_with_cpptraj(pdb_path)
        num_close_atoms = calculate_close_atoms_pymol(pdb_path)

        data.append([os.path.basename(pdb_path), hbond_between, hbond_combined, hbond_cap, hbond_leeway, num_close_atoms])

        # Check criteria and copy file if it meets the conditions
        if hbond_between >= MIN_HBONDS and num_close_atoms >= MIN_CLOSE_ATOMS:
            shutil.copy(pdb_path, f"{filter_folder}")

    # Create DataFrame and save results
    df = pd.DataFrame(data, columns=["PDB_File", "H-Bonds_Between", "H-Bonds_Combined", "H-Bonds_Cap", "H-Bonds_Leeway", "Close_Atoms"])
    df.to_csv(results_file, index=False, sep='\t')
    return df

# 8. Plotting functions
def plot_histograms(df):
    plt.figure(figsize=(10, 10))
    n, bins, patches = plt.hist(df["H-Bonds_Between"], bins=10, alpha=0.7, color='orange')
    plt.xlabel("HB") # Between cap and leeway
    plt.ylabel("Frequency")
    plt.title("HB between Cap and Leeway")
    for i in range(len(n)):
        plt.text(bins[i] + (bins[i+1] - bins[i]) / 2, n[i], str(int(n[i])), 
                 ha='center', va='bottom', fontsize=25)
    plt.tight_layout() 
    plt.savefig(hbond_histogram_file)
    plt.show()

    plt.figure(figsize=(10, 10))
    n, bins, patches = plt.hist(df["Close_Atoms"], bins=10, alpha=0.7, color='orange')
    plt.xlabel("Num of atoms")
    plt.ylabel("Frequency")
    plt.title("Number of 3Å Atoms")
    for i in range(len(n)):
        plt.text(bins[i] + (bins[i+1] - bins[i]) / 2, n[i], str(int(n[i])), 
                 ha='center', va='bottom', fontsize=25)
    plt.tight_layout()
    plt.savefig(close_atom_histogram_file)
    plt.show()

def plot_scatter(df):
    plt.figure(figsize=(10, 10))
    plt.scatter(df["H-Bonds_Between"], df["Close_Atoms"], alpha=1, color="blue", s=200)
    plt.xlabel("HB between cap and leeway") # Between cap and leeway
    plt.ylabel("Num of atoms")
    plt.title("HB vs. Atoms in 3Å")
    plt.tight_layout()
    plt.savefig(scatter_plot_file)
    plt.show()

# Main script to run everything
if __name__ == "__main__":
    
    # Get list of PDB files, converting CIF files if needed
    pdb_files = get_pdb_files(input_folder)
    pdb_files = add_hydrogens(pdb_files)

    # Calculate hydrogen bonds and close atoms, then filter and save results
    df = filter_and_copy_pdbs(pdb_files)

    # Generate histograms and scatter plot
    plot_histograms(df)
    plot_scatter(df)
