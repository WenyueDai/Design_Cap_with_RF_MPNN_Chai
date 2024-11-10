#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH --gres=gpu:rtx2080:1
#SBATCH -c 2
#SBATCH --output=example_1.out

# Activate environment (use conda or adjust as needed)
source activate mlfold

# Input and output directories
folder_with_pdbs="/home/eva/0_bury_charged_pair/0_With_initial_pdb/1_interesting_start/3ult_HDHD"
output_dir="/home/eva/0_bury_charged_pair/0_With_initial_pdb/2_MPNN/3ult"
mkdir -p $output_dir

# Paths for intermediate files
path_for_parsed_chains="$output_dir/parsed_pdbs.jsonl"
path_for_assigned_chains="$output_dir/assigned_pdbs.jsonl"
path_for_fixed_positions="$output_dir/fixed_pdbs.jsonl"
chains_to_design="A"
fixed_residue_indices="43 58 72 86"  # Replace with the specific indices you want to keep fixed

# Step 1: Parse PDB files to extract sequences and residue indices
python /home/eva/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

# Step 3: Assign chains for design (i.e., which chains will be designed)
python /home/eva/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

# Step 4: Create a fixed positions dictionary using the provided indices
python /home/eva/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_residue_indices"

# Step 5: Run ProteinMPNN with fixed positions
python /home/eva/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 500 \
        --sampling_temp "0.1" \
        --omit_AAs "CWYF" \
        --seed 37 \
        --batch_size 1

echo "ProteinMPNN run complete with specified residues fixed."

