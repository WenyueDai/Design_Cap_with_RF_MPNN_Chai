import os
import pandas as pd
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqUtils import seq1
from MDAnalysis import Universe
import warnings

"""After RFdiffusion, there are some unsatisfied cap that has little chance of being a good cap after MPNN.
These cap include the ones that are far away from the main body, have a severe angle and has a too close RMSD 
with the main body because that would make it more like an extension rather than a cap. I filter them out and 
only keep the good one for MPNN.

It output either Ncap-mainbody or Ccap-mainbody

As long as the RMSD, angle and distance histogram and 2D plot was generated from last script.
It is still not very clear what is the 'good cap limitation' for all three parameters.
What I do is to first slice based on angle (0,5) (5, 10), (10, 15), (15, 20), (20, 25) and then run the script
to seperate out the pdb based on the sliced. Then I manually look at the pdb structure filtered out in each slices.
Ater find out the good range, use it to slice distance (18, 22), (22, 23), (23, 24), (24, 25), (25>)
Then find out the good range, use both parameter to slice rmsd (5, 8), (8, 9), (9,10), (10, 11), (11,15)
After applying the narrowed range, I in total filter out around 200 pdb structure for each N cap and C cap for m7 and m8
"""

# Function to filter PDB files based on given ranges for N and C cap parameters
def filter_pdb_files(parameter_file, n_rmsd_range, n_distance_range, n_angle_range,
                     c_rmsd_range, c_distance_range, c_angle_range, output_folder):
    # Load the parameter file into a DataFrame
    df = pd.read_csv(parameter_file, sep='\t')
    
    # Create subfolders for N cap and C cap within the output folder
    n_cap_folder = os.path.join(output_folder, "N_cap")
    c_cap_folder = os.path.join(output_folder, "C_cap")
    os.makedirs(n_cap_folder, exist_ok=True)
    os.makedirs(c_cap_folder, exist_ok=True)
    
    # Filter for PDB files that meet the conditions for N cap and C cap
    filtered_n_cap = df[(df['N_Cap_RMSD'] >= n_rmsd_range[0]) & (df['N_Cap_RMSD'] <= n_rmsd_range[1]) &
                        (df['N_Cap_Distance'] >= n_distance_range[0]) & (df['N_Cap_Distance'] <= n_distance_range[1]) &
                        (df['N_Cap_Angle'] >= n_angle_range[0]) & (df['N_Cap_Angle'] <= n_angle_range[1])]
    
    filtered_c_cap = df[(df['C_Cap_RMSD'] >= c_rmsd_range[0]) & (df['C_Cap_RMSD'] <= c_rmsd_range[1]) &
                        (df['C_Cap_Distance'] >= c_distance_range[0]) & (df['C_Cap_Distance'] <= c_distance_range[1]) &
                        (df['C_Cap_Angle'] >= c_angle_range[0]) & (df['C_Cap_Angle'] <= c_angle_range[1])]
    
    print(f"Number of PDB files that meet N cap criteria: {len(filtered_n_cap)}")
    print(f"Number of PDB files that meet C cap criteria: {len(filtered_c_cap)}")
    
    # Process and copy files based on N cap filtering
    for _, row in filtered_n_cap.iterrows():
        pdb_path = row['PDB_File']
        if os.path.isfile(pdb_path):
            print(f"Processing N cap for PDB file: {pdb_path}")
            process_and_save_filtered_structure(pdb_path, "N", n_cap_folder)
    
    # Process and copy files based on C cap filtering
    for _, row in filtered_c_cap.iterrows():
        pdb_path = row['PDB_File']
        if os.path.isfile(pdb_path):
            print(f"Processing C cap for PDB file: {pdb_path}")
            process_and_save_filtered_structure(pdb_path, "C", c_cap_folder)

# Function to process PDB structure and retain only the analyzed cap and main body
def process_and_save_filtered_structure(pdb_path, cap_type, output_folder):
    try:
        # Extract sequence to determine cap lengths
        sequence = extract_sequence_from_pdb(pdb_path)
        print(f"Extracted sequence from {pdb_path}: {sequence}")

        # Determine N cap and C cap lengths using glycine detection
        n_cap_len, c_cap_len = determine_glycine_caps(sequence)
        print(f"Determined N cap length: {n_cap_len}, C cap length: {c_cap_len} for {pdb_path}")

        if n_cap_len == 0 or c_cap_len == 0:
            print(f"Skipping {pdb_path} due to insufficient N or C cap length.")
            return

        # Load PDB file with MDAnalysis
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = Universe(pdb_path)

        # Define atom groups for N cap, C cap, and main body
        n_cap = u.select_atoms(f"resid 1:{n_cap_len}")
        c_cap = u.select_atoms(f"resid {len(sequence) - c_cap_len + 1}:{len(sequence)}")
        main_body = u.select_atoms(f"resid {n_cap_len + 1}:{len(sequence) - c_cap_len}")

        # Debug information for atom group selection
        print(f"N cap atom count for {pdb_path}: {len(n_cap)}")
        print(f"C cap atom count for {pdb_path}: {len(c_cap)}")
        print(f"Main body atom count for {pdb_path}: {len(main_body)}")

        # Combine cap with main body depending on the cap type
        if cap_type == "N":
            cap_and_main_body = n_cap + main_body
            output_file = os.path.join(output_folder, f"Ncap_{os.path.basename(pdb_path)}")
        elif cap_type == "C":
            cap_and_main_body = c_cap + main_body
            output_file = os.path.join(output_folder, f"Ccap_{os.path.basename(pdb_path)}")

        # Debug: Print the total number of atoms to be written
        print(f"Total atom count for {cap_type} cap + main body in {pdb_path}: {len(cap_and_main_body)}")

        # Save the filtered structure with only the analyzed cap and main body
        if len(cap_and_main_body) > 0:
            cap_and_main_body.write(output_file)
            print(f"Filtered PDB saved to {output_file}")
        else:
            print(f"No atoms found to save for {cap_type} cap in {pdb_path}")

    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")

# Function to extract protein sequence from PDB using Biopython
def extract_sequence_from_pdb(pdb_path):
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_path)
        ppb = PPBuilder()
        sequence = ""
        for pp in ppb.build_peptides(structure):
            sequence += str(pp.get_sequence())
        return sequence
    except Exception as e:
        print(f"Error extracting sequence from {pdb_path}: {e}")
        return ""

# Helper function to determine glycine-based N cap and C cap lengths
def determine_glycine_caps(sequence):
    n_cap_len, c_cap_len = 0, 0
    # Detect N cap (continuous glycine residues at the start)
    for i, res in enumerate(sequence):
        if res == 'G':
            n_cap_len += 1
        else:
            break
    # Detect C cap (continuous glycine residues at the end)
    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] == 'G':
            c_cap_len += 1
        else:
            break
    return n_cap_len, c_cap_len


# Example usage
parameter_file = "/home/eva/0_bury_charged_pair/5_Pipeline/20241101_rfdiffusion_m7/1_cap_calculation/cap_parameters.txt"
n_rmsd_range = (5.0, 11.0)          # Min and Max RMSD for N cap
n_distance_range = (22.0, 40.0)     # Min and Max Distance for N cap
n_angle_range = (0.0, 20.0)         # Min and Max Angle for N cap
c_rmsd_range = (5.0, 11.0)          # Min and Max RMSD for C cap
c_distance_range = (23.0, 26.0)     # Min and Max Distance for C cap
c_angle_range = (0.0, 15.0)         # Min and Max Angle for C cap
output_folder = "/home/eva/0_bury_charged_pair/5_Pipeline/20241101_rfdiffusion_m7/1_cap_calculation"

# Ensure the output folder and subfolders exist
os.makedirs(output_folder, exist_ok=True)

# Run the filtering function
filter_pdb_files(parameter_file, n_rmsd_range, n_distance_range, n_angle_range,
                 c_rmsd_range, c_distance_range, c_angle_range, output_folder)
