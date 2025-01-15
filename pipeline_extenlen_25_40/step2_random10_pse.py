import os
import random
import pymol
from pymol import cmd

"""This is the script I use after step2_filter_s2_move.py, to randomly select 10 pdb structure and generate pse file to look for a give pdb folder path
"""

def select_and_save_pse(folder_path, pse_name="combined_pdbs.pse"):
    # Check if the provided folder path exists
    if not os.path.isdir(folder_path):
        print(f"The folder path '{folder_path}' does not exist.")
        return

    # List all PDB files in the folder
    pdb_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
    if len(pdb_files) < 10:
        print("Not enough PDB files in the folder to select 10. Please add more PDB files.")
        return

    # Randomly select 10 PDB files
    selected_pdb_files = random.sample(pdb_files, 10)
    print("Selected PDB files:", selected_pdb_files)

    # Start PyMOL in command-line mode
    pymol.finish_launching(["pymol", "-cq"])

    # Load each selected PDB file into PyMOL
    for pdb_file in selected_pdb_files:
        pdb_path = os.path.join(folder_path, pdb_file)
        cmd.load(pdb_path)

    # Save the session as a .pse file
    cmd.save(os.path.join(folder_path, pse_name))
    print(f"Saved PyMOL session as '{pse_name}' in '{folder_path}'.")

    # Quit PyMOL session
    cmd.quit()

# Example usage:
folder_path = "/home/eva/0_bury_charged_pair/5_Pipeline/20241101_rfdiffusion_m8/1_cap_calculation/11_15_N_cap"  # Replace with your folder path
select_and_save_pse(folder_path)
