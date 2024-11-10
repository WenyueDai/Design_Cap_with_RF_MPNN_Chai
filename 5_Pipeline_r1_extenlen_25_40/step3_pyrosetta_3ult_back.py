import os
import pyrosetta

pyrosetta.init("--ex1 --ex2")

"""After filter out the good protein scaffold, I use pyrosetta to map the 3ult sequence back.
That is just need to find out the region that is not continuous glycine - that should be
the main body made of valine and glutamamte.
"""

def one_letter_to_three(letter):
    one_to_three = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }
    return one_to_three.get(letter)

def detect_glycine_caps(pose):
    """Automatically detect the N cap and C cap lengths based on glycine residues in PyRosetta."""
    n_cap_len, c_cap_len = 0, 0
    # Detect N cap (continuous glycine residues at the start)
    for i in range(1, pose.total_residue() + 1):
        if pose.residue(i).name1() == 'G':
            n_cap_len += 1
        else:
            break
    # Detect C cap (continuous glycine residues at the end)
    for i in range(pose.total_residue(), 0, -1):
        if pose.residue(i).name1() == 'G':
            c_cap_len += 1
        else:
            break
    return n_cap_len, c_cap_len

def modify_main_body_sequence(pose, input_sequence, n_cap_len, c_cap_len):
    """Replace only the main body sequence, leaving the N and C caps unchanged."""
    main_body_start = n_cap_len + 1
    main_body_end = pose.total_residue() - c_cap_len

    # Calculate the expected main body length
    expected_main_body_length = main_body_end - main_body_start + 1
    print(f"Detected N cap length: {n_cap_len}, C cap length: {c_cap_len}")
    print(f"Main body start: {main_body_start}, Main body end: {main_body_end}")
    print(f"Expected main body length: {expected_main_body_length}")
    print(f"Input sequence length: {len(input_sequence)}")

    # Check if input sequence length matches the main body length
    if len(input_sequence) != expected_main_body_length:
        raise ValueError("Input sequence length does not match the number of residues in the main body.")

    # Replace residues in the main body
    for i, resi in enumerate(range(main_body_start, main_body_end + 1)):
        target_residue = input_sequence[i]
        target_residue_3letter = one_letter_to_three(target_residue)

        # Only mutate if the target residue is different from the current residue
        if pose.residue(resi).name1() != target_residue:
            pose.replace_residue(
                resi, 
                pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(
                    pose.residue_type_set_for_pose().name_map(target_residue_3letter)
                ), 
                True
            )

def process_folder(input_folder, output_folder, input_sequence):
    """Process each PDB file in the input folder, modify the main body sequence, and save it to output folder."""
    os.makedirs(output_folder, exist_ok=True)

    for pdb_file in os.listdir(input_folder):
        if pdb_file.endswith(".pdb"):
            input_pdb_path = os.path.join(input_folder, pdb_file)
            output_pdb_path = os.path.join(output_folder, f"modified_{pdb_file}")

            # Load the pose and detect cap regions
            pose = pyrosetta.pose_from_file(input_pdb_path)
            n_cap_len, c_cap_len = detect_glycine_caps(pose)

            # Modify only the main body sequence
            modify_main_body_sequence(pose, input_sequence, n_cap_len, c_cap_len)

            # Save the modified pose
            pose.dump_pdb(output_pdb_path)
            print(f"Processed {pdb_file} -> {output_pdb_path}")

# Example usage
n_cap_folder = '/home/eva/0_bury_charged_pair/5_Pipeline/20241101_rfdiffusion_m7/1_cap_calculation/N_cap'
c_cap_folder = '/home/eva/0_bury_charged_pair/5_Pipeline/20241101_rfdiffusion_m7/1_cap_calculation/C_cap'
output_n_folder = '/home/eva/0_bury_charged_pair/5_Pipeline/20241101_rfdiffusion_m7/1_cap_calculation/N_cap_add_3ultback'
output_c_folder = '/home/eva/0_bury_charged_pair/5_Pipeline/20241101_rfdiffusion_m7/1_cap_calculation/C_cap_add_3ultback'

# Define the sequence to replace the main body with
main_body_sequence = "PNTISGSNNTVRSGSKNVLAGNDNTVISGDNNSVSGSNNTVVSGNDNTVTGSNHVVSGTNHIVTDNNNNVSGNDNNVSGSFHTVSGGHNTVSGSNNTVSGSNHVVSGSNKVVTD"  # Replace with the actual sequence

# Process N cap and C cap folders
process_folder(n_cap_folder, output_n_folder, main_body_sequence)
process_folder(c_cap_folder, output_c_folder, main_body_sequence)
