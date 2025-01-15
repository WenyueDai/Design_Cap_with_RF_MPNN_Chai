#!/bin/bash

# Input: Folder containing FASTA files
input_folder="/home/eva/0_bury_charged_pair/5_Pipeline/20241101_rfdiffusion_m7/1_cap_calculation/N_cap_add_3ultback/mpnn/seqs"

# Check if folder was provided
if [ -z "$input_folder" ]; then
  echo "Usage: $0 <input_folder>"
  exit 1
fi

# Create output directory if it doesn't exist
output_dir="/home/eva/0_bury_charged_pair/5_Pipeline/20241103_mpnn_m7/N_cap"
mkdir -p "$output_dir"

# Loop through each FASTA file in the input folder
for parent_fasta in "$input_folder"/*.fa; do

  # Extract the base name of the parent FASTA (without the extension)
  base_name=$(basename "$parent_fasta" .fasta)

  # Initialize sample number counter
  sample_num=0

  # Read the parent FASTA file line by line
  while IFS= read -r line; do
    # Check if the line is a header (starts with '>')
    if [[ $line == ">"* ]]; then
      # Extract sample number and global score from the header line
      sample_num=$(echo "$line" | grep -oP 'sample=\K\d+' || echo "$sample_num")
      global_score=$(echo "$line" | grep -oP 'global_score=\K[\d\.]+' || echo "unknown")

      # Prepare the new header in the required format
      new_header=">sample=${sample_num}_global_score=${global_score}"

      # Create a new FASTA file for this sample
      output_file="${output_dir}/${base_name}_sample${sample_num}.fasta"
      echo "$new_header" > "$output_file"
    else
      # If it's a sequence line, append it to the current output file
      echo "$line" >> "$output_file"
    fi
  done < "$parent_fasta"

done
# Remove files containing "sample0" in the name from the output directory
find "$output_dir" -type f -name "*sample0*.fasta" -delete
echo "FASTA extraction complete. Files saved in $output_dir/"
