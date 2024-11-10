import os
import re
import subprocess
from Bio import AlignIO
from math import log2

# Paths to CD-HIT and MAFFT executables
CD_HIT_EXECUTABLE = "cd-hit"  # Ensure CD-HIT is installed and accessible
MAFFT_EXECUTABLE = "mafft"    # Ensure MAFFT is installed and accessible

# Define input and output folder paths
input_fasta = "/home/eva/0_bury_charged_pair/0_With_initial_pdb/2_MPNN/3ult/seqs/3ult_HDHD.fa"
cleaned_fasta = "/home/eva/0_bury_charged_pair/0_With_initial_pdb/2_MPNN/3ult/seqs/3ult_HDHD_clean.fa"
aligned_fasta = "/home/eva/0_bury_charged_pair/0_With_initial_pdb/2_MPNN/3ult/seqs/aligned_sequences.fasta"
output_folder = "/home/eva/0_bury_charged_pair/0_With_initial_pdb/2_MPNN/3ult/top100"
conservation_scores = os.path.join(output_folder, "conservation_scores.txt")

# Function to clean headers in the FASTA file
def clean_fasta_headers(input_fasta, output_fasta):
    """
    Clean the headers in the input FASTA file by extracting sample number and global score.
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith(">"):
                # Extract sample number and global score
                sample_num = re.search(r'sample=(\d+)', line)
                global_score = re.search(r'global_score=([\d\.]+)', line)
                sample_num = sample_num.group(1) if sample_num else "unknown"
                global_score = global_score.group(1) if global_score else "unknown"
                new_header = f">sample={sample_num}_global_score={global_score}"
                outfile.write(new_header + "\n")
            else:
                outfile.write(line)
    print(f"FASTA headers cleaned and saved to: {output_fasta}")

# Function to filter sequences using CD-HIT (for diversity)
def filter_sequences_with_cd_hit(input_fasta, output_fasta, identity_threshold=0.90):
    """
    Filter out repetitive sequences using CD-HIT on a cleaned FASTA file.
    """
    try:
        print(f"Running CD-HIT on {input_fasta}...")
        subprocess.run([
            CD_HIT_EXECUTABLE,
            "-i", input_fasta,
            "-o", output_fasta,
            "-c", str(identity_threshold),
            "-n", "5",  # Use word size 5 for faster processing
            "-G", "0",
            "-aL", "0.5",
            "-d", "0"
        ], check=True)
        print(f"Filtered unique sequences saved to: {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {input_fasta}: {e}")

# Function to align sequences using MAFFT with high accuracy
def align_sequences(input_fasta, aligned_fasta, threads=4):
    """
    Align sequences using MAFFT with high accuracy and save to an output file.
    """
    try:
        print(f"Aligning sequences from {input_fasta} using MAFFT with high accuracy...")
        subprocess.run([
            MAFFT_EXECUTABLE,
            "--maxiterate", "1000",
            "--localpair",
            "--thread", str(threads),
            input_fasta
        ], stdout=open(aligned_fasta, 'w'), stderr=subprocess.DEVNULL)
        print(f"Aligned sequences saved to: {aligned_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error aligning sequences: {e}")

# Function to calculate sequence entropy for conservation
def calculate_sequence_entropy(alignment):
    """
    Calculate sequence entropy for each column of an alignment.
    """
    entropy_scores = []
    for col in range(len(alignment[0])):
        column = [str(record.seq[col]) for record in alignment if record.seq[col] != "-"]
        if len(column) == 0:
            entropy_scores.append(0)
            continue
        freq = {base: column.count(base) / len(column) for base in set(column)}
        entropy = -sum(f * log2(f) for f in freq.values() if f > 0)
        entropy_scores.append(entropy)
    return entropy_scores

# Function to calculate conservation scores using entropy
def calculate_conservation_scores(fasta_file, output_conservation_scores):
    """
    Calculate conservation scores for the sequences in the provided aligned FASTA file using entropy.
    """
    alignment = AlignIO.read(fasta_file, "fasta")
    entropy_scores = calculate_sequence_entropy(alignment)
    with open(output_conservation_scores, 'w') as outfile:
        for i, score in enumerate(entropy_scores):
            outfile.write(f"Position {i+1}: {score:.4f}\n")
    print(f"Conservation scores saved to: {output_conservation_scores}")

# Function to select top sequences based on diversity and conservation
def select_top_sequences_by_diversity_and_conservation(output_fasta, n=100, output_folder="selected_sequences", combined_output_fasta="top_100_sequences.fasta"):
    """
    Select top sequences based on diversity (from CD-HIT) and conservation scores.
    """
    # Read sequences from the CD-HIT output
    sequences = []
    with open(output_fasta, 'r') as f:
        seq = ''
        header = ''
        for line in f:
            if line.startswith('>'):
                if seq != '':
                    sequences.append((header, seq))
                    seq = ''
                header = line.strip()
            else:
                seq += line.strip()
        if seq != '':
            sequences.append((header, seq))

    # If there are fewer sequences than n, adjust n
    n = min(n, len(sequences))

    # Select top N sequences
    selected_sequences = sequences[:n]

    # Write combined and individual FASTA files
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    combined_output_fasta = os.path.join(output_folder, combined_output_fasta)

    with open(combined_output_fasta, 'w') as combined_outfile:
        for idx, (header, seq) in enumerate(selected_sequences):
            combined_outfile.write(f"{header}\n{seq}\n")
            individual_fasta = os.path.join(output_folder, f"sequence_{idx+1}.fasta")
            with open(individual_fasta, 'w') as individual_outfile:
                individual_outfile.write(f"{header}\n{seq}\n")
    print(f"Top {n} sequences saved to {combined_output_fasta} and {output_folder}")

# Main processing block
if __name__ == "__main__":
    # Ensure output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Step 1: Clean FASTA headers
    clean_fasta_headers(input_fasta, cleaned_fasta)

    # Step 2: Run CD-HIT on the cleaned FASTA file (for diversity)
    output_fasta = os.path.join(output_folder, "unique_combined_sequences.fasta")
    filter_sequences_with_cd_hit(cleaned_fasta, output_fasta)

    # Step 3: Align the sequences using MAFFT
    align_sequences(output_fasta, aligned_fasta)

    # Step 4: Calculate conservation scores using entropy
    calculate_conservation_scores(aligned_fasta, conservation_scores)

    # Step 5: Select top sequences based on diversity
    select_top_sequences_by_diversity_and_conservation(output_fasta, n=100, output_folder=output_folder)

    print(f"\nProcessing complete. Filtered sequences can be found in the specified output folder.")
