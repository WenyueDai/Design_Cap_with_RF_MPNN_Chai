import pyrosetta # type: ignore
from pymol import cmd # type: ignore
pyrosetta.init("--ex1 --ex2")

"""This script perform the following steps:
1. Clean up 3ult.pdb and keep only chain A.
2. Repack 3ult.pdb with a combination of 10*V + 94*E + 10*V as a starting structure for RFdiffusion. 
"""

#conda activate getcontact

def clean_and_renumber_pdb(input_file, output_file, chains_to_keep="A"):
    """
    Cleans and renumbers a PDB file:
    - Removes water and common ions (NA, CL, MG, CA, ZN).
    - Keeps only atoms with alternate locations 'A' or ' '.
    - Filters to keep only specified chain(s), default is chain A.
    - Renumbers residues and atoms sequentially
    """
    cmd.load(input_file, "protein")
    cmd.remove("resn HOH or resn NA+CL+MG+CA+ZN or hetatm")
    cmd.remove("not alt ''+A")
    if chains_to_keep:
        cmd.remove(f"not chain {chains_to_keep}")

    # Renumber residues to start from 1
    starting_residue = int(cmd.get_model("all").atom[0].resi)  # Find the starting residue number
    cmd.alter("all", f"resi=str(int(resi) - {starting_residue - 1})")  # Adjust residues to start from 1

    # Renumber atoms sequentially starting from 1
    atom_counter = 1
    cmd.alter("all", "serial=atom_counter; atom_counter+=1", space={'atom_counter': atom_counter})
    cmd.sort()
    cmd.save(output_file, "protein")
    cmd.delete("all")
    return "PDB clean up finished"

def parse_pdb_sequence(file_path, chain_id='A'):
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    sequence = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(('ATOM')) and line[21].strip() == chain_id:
                res_name = line[17:20].strip()
                if res_name in three_to_one:
                    sequence[line[22:27].strip()] = three_to_one[res_name]   
    return "".join(sequence.values())

#Sequence of chain A: PNTISGSNTVRSGSKNVLAGNDNTVISGDNSVSGSNTVSGNDNTVTGSNHVSGTNHIVTDNVSGNDNVSGSFHTVSGHNTVSGSNTVSGSNHVSGSNKVTD
# Provide the path to your PDB file
pdb_file_path = '/home/eva/0_bury_charged_pair/2_No_initial_de_novo/2_rfdiffusion_cap/input/3ult_cleaned_renumbered_monomer.pdb'

sequence = parse_pdb_sequence(pdb_file_path)
print(f"Sequence of chain A: {sequence}")

def one_letter_to_three(letter):
    one_to_three = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }
    return one_to_three.get(letter)

def modify_protein_sequence(input_pdb: str, output_pdb: str, input_sequence: str) -> None:
    # Initialize PyRosetta
    #pyrosetta.init("--ex1 --ex2")
    
    # Load the pose from the PDB file
    pose = pyrosetta.pose_from_file(input_pdb)
    
    # Check if the input sequence length matches the total residues in the PDB file
    if len(input_sequence) != pose.total_residue():
        raise ValueError("Input sequence length does not match the number of residues in the PDB file.")
    
    # Map the input sequence to the PDB structure
    for resi in range(1, pose.total_residue() + 1):
        # Get the target residue type from the input sequence
        target_residue = input_sequence[resi - 1]
        # Convert the one-letter code to three-letter code
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
    
    # Repack the structure
    scorefxn = pyrosetta.rosetta.core.scoring.get_score_function()
    pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
    pack_mover.task_factory(tf)
    pack_mover.apply(pose)
    
    # Dump the modified pose to the defined file
    pose.dump_pdb(output_pdb)

# Example usage
start_pdb = "/home/eva/0_bury_charged_pair/5_Pipeline_r1_extenlen_25_40/input_pdb/3ult.pdb"
dammy_seq = 'EVE'*1 + "V"*9 + 'E'*90 + 'V'*8 + 'EVEE'*1
input_pdb = '/home/eva/0_bury_charged_pair/5_Pipeline_r1_extenlen_25_40/input_pdb/3ult_cleaned_renumbered_monomer.pdb'
output_pdb = '/home/eva/0_bury_charged_pair/5_Pipeline_r1_extenlen_25_40/input_pdb/3ult_cleaned_renumbered_monomer_VE.pdb'
input_sequence = dammy_seq  # Replace with the actual input sequence
clean_and_renumber_pdb(start_pdb, input_pdb, chains_to_keep="A")
modify_protein_sequence(input_pdb, output_pdb, input_sequence)
