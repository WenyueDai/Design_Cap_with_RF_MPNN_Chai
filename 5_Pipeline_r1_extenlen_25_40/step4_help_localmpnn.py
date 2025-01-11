import os
import subprocess
import json

"""This is solublempnn to find the sequence for both capping, and the 14AA extended from the capping sequence
"""

#conda activate mlfold

# Input and output directories
folder_with_pdbs = "/home/eva/0_bury_charged_pair/5_Pipeline/20241101_rfdiffusion_m7/1_cap_calculation/N_cap_add_3ultback"
output_dir = "/home/eva/0_bury_charged_pair/5_Pipeline/20241101_rfdiffusion_m7/1_cap_calculation/N_cap_add_3ultback/mpnn"

# Ensure the output folder exists
os.makedirs(output_dir, exist_ok=True)

# Paths for intermediate files
path_for_parsed_chains = f"{output_dir}/parsed_pdbs.jsonl"
path_for_assigned_chains = f"{output_dir}/assigned_pdbs.jsonl"
path_for_fixed_positions = f"{output_dir}/fixed_positions_caps.jsonl"
chains_to_design = "A"

# Step 1: Parse PDB files to extract sequences and residue indices
# Command to call parse_multiple_chains.py script
command_parse_chains = [
    "python",
    "/home/eva/ProteinMPNN/helper_scripts/parse_multiple_chains.py",
    "--input_path", folder_with_pdbs,
    "--output_path", path_for_parsed_chains
]

# Run the parsing command
print("Running parse_multiple_chains.py to parse PDB chains...")
subprocess.run(command_parse_chains, check=True)
print("Parsing complete.")

# Step 2: Identify N and C Cap regions and set fixed positions for main body
cap_positions_dict = {}
assigned_chains_dict = {}

# Helper to detect N and C caps based on continuous glycine residues
def identify_cap_regions(sequence):
    n_cap_end_resi = 0
    c_cap_start_resi = len(sequence)

    # Detect N cap: continuous glycine from the start
    for res in sequence:
        if res == 'G':
            n_cap_end_resi += 1
        else:
            break

    # Detect C cap: continuous glycine from the end
    for res in sequence[::-1]:
        if res == 'G':
            c_cap_start_resi -= 1
        else:
            break

    leeway = 14
    # Determine if only N cap or C cap is present
    if n_cap_end_resi > 0 and c_cap_start_resi == len(sequence):
        return list(range(1, n_cap_end_resi + leeway + 1)), []  # Only N cap, using residue indexes
    elif c_cap_start_resi < len(sequence) and n_cap_end_resi == 0:
        return [], list(range(c_cap_start_resi - leeway, len(sequence) + 1))  # Only C cap, using residue indexes
    else:
        return [], []

# Read parsed JSON and extract positions to fix main body
with open(path_for_parsed_chains, 'r') as json_file:
    parsed_data = [json.loads(line) for line in json_file]

# Generate fixed positions for each chain
for entry in parsed_data:
    pdb_name = entry["name"].split('.')[0]  # Remove file extension if present
    sequence_key = f"seq_chain_{chains_to_design}"

    # Ensure the chain is present in the parsed entry
    if sequence_key not in entry:
        print(f"Warning: Chain {chains_to_design} not found in entry {pdb_name}, skipping.")
        continue

    sequence = entry[sequence_key]

    # Identify cap regions
    n_cap, c_cap = identify_cap_regions(sequence)

    # Define main body fixed positions
    main_body_fixed_positions = [i for i in range(1, len(sequence)+1) if i not in (n_cap+c_cap)]
    
    # Store in dictionary format required by ProteinMPNN
    cap_positions_dict[pdb_name] = {chains_to_design: main_body_fixed_positions}
    print(f"Fixed positions for {pdb_name}: {cap_positions_dict[pdb_name]}")

    # Add assigned chain information with two lists as expected by ProteinMPNN
    assigned_chains_dict[pdb_name] = [[chains_to_design], []]  # Masked chains, visible chains

with open(f"{output_dir}/fixed_positions_caps.jsonl", 'w') as file:
    file.write(json.dumps(cap_positions_dict))
    
# Write assigned chains to JSONL format using the helper script
with open(f"{output_dir}/temp_assigned_chains.json", 'w') as temp_file:
    for pdb, chains in assigned_chains_dict.items():
        temp_file.write(json.dumps({pdb: chains}) + '\n')

command_assign_chains = [
    "python",
    "/home/eva/ProteinMPNN/helper_scripts/assign_fixed_chains.py",
    "--input_path", path_for_parsed_chains,
    "--output_path", path_for_assigned_chains,
    "--chain_list", chains_to_design
]
subprocess.run(command_assign_chains, check=True)
print("Assigned chains JSONL generated.")

print("Generated assigned chains and fixed positions JSON files.")

# Step 3: Run ProteinMPNN with generated JSON files
# Command to call protein_mpnn_run.py script
command_mpnn_run = [
    "python",
    "/home/eva/ProteinMPNN/protein_mpnn_run.py",
    "--jsonl_path", path_for_parsed_chains,
    "--chain_id_jsonl", path_for_assigned_chains,
    "--fixed_positions_jsonl", path_for_fixed_positions,
    "--out_folder", output_dir,
    "--use_soluble_model",
    "--num_seq_per_target", "2",
    "--sampling_temp", "0.1",
    "--omit_AAs", "CWY",
    "--seed", "37",
    "--batch_size", "1"
]

# Run the ProteinMPNN command
print("Running protein_mpnn_run.py to design sequences...")
subprocess.run(command_mpnn_run, check=True)
print("ProteinMPNN run complete.")
