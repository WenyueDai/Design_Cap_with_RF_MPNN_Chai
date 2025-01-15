import os
import Bio.PDB
import pymol
from pymol import cmd, stored
import shutil

# Define the folders for input and output PDB files
pdb_folder = "/home/eva/0_RFdiffusion_MPNN_alphafold_pip_ulti/20240912_keep_2D/A1_native_pdb_clean"
cleaned_pdb_folder = "/home/eva/0_bury_charged_pair/1_No_initial_pdb/1_find_initial_pdb/1_cleaned_solenoid_pdb"
buried_charged_pairs_folder = "/home/eva/0_bury_charged_pair/1_No_initial_pdb/1_find_initial_pdb/2_buried_charge_solenoid_pdb"
os.makedirs(cleaned_pdb_folder, exist_ok=True)
os.makedirs(buried_charged_pairs_folder, exist_ok=True)

# Initialize PyMOL in quiet mode
pymol.finish_launching(['pymol', '-cq'])

# Define a class to selectively save only protein atoms (remove water and metal ions)
class CleanStructure(Bio.PDB.Select):
    def accept_residue(self, residue):
        if residue.id[0] != ' ' or residue.get_resname() in ['HOH', 'ZN', 'MG', 'CA', 'NA', 'FE']:
            return 0  # Reject water and metal ions
        return 1  # Accept all other residues

# Function to clean the PDB file and remove water and metal ions
def clean_pdb(pdb_file, cleaned_pdb_file):
    with open(pdb_file, 'r') as infile, open(cleaned_pdb_file, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM"):
                break
            outfile.write(line)
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)
        io = Bio.PDB.PDBIO()
        io.set_structure(structure)
        io.save(outfile, select=CleanStructure())  
    print(f"Cleaned PDB saved to {cleaned_pdb_file}")

# Function to calculate SASA for charged residues and check for oppositely charged pairs
def calculate_sasa_for_charged_pairs(pdb_file):
    print(f"[INFO] Loading {pdb_file} in PyMOL...")
    cmd.load(pdb_file)
    
    # Set PyMOL parameters for SASA calculation
    cmd.set('dot_solvent', 1)
    cmd.set('dot_density', 3)

    # Store unique charged residues in a dictionary
    stored.charged_residues = {}
    cmd.iterate('resn ASP+GLU+LYS+ARG+HIS', 'stored.charged_residues[(chain, resi)] = resn')
    
    # List of opposite charged residue pairs
    opposites = {'ASP': ['LYS', 'ARG', 'HIS'], 'GLU': ['LYS', 'ARG', 'HIS'], 
                 'LYS': ['ASP', 'GLU'], 'ARG': ['ASP', 'GLU'], 'HIS': ['ASP', 'GLU']}
    
    # Search for low SASA and opposite charged pairs
    charged_pairs = []
    for (chain1, resi1), resn1 in stored.charged_residues.items():
        sasa1 = cmd.get_area(f"resi {resi1} and chain {chain1} and not name N+CA+C+O")
        if sasa1 < 1.0:
            for (chain2, resi2), resn2 in stored.charged_residues.items():
                if (resn1 in opposites) and (resn2 in opposites[resn1]) and sasa1 < 1.0:
                    distance = cmd.distance(f"dist_{resi1}_{resi2}", f"resi {resi1} and chain {chain1}",
                                            f"resi {resi2} and chain {chain2}")
                    if distance < 10.0:
                        print(f"[INFO] Found buried opposite charged pair: {resn1} {resi1}-{resn2} {resi2}")
                        charged_pairs.append((resi1, resn1, resi2, resn2))
    
    cmd.delete('all')  # Clean up PyMOL for the next PDB file
    return charged_pairs

# Function to copy PDB files with buried opposite charged pairs
def copy_pdb_with_charged_pairs(pdb_file, charged_pairs, output_folder):
    if charged_pairs:
        destination = os.path.join(output_folder, os.path.basename(pdb_file))
        shutil.copy(pdb_file, destination)
        print(f"[INFO] Copied {pdb_file} to {output_folder}")

# Function to generate PyMOL script to highlight charged pairs
def generate_pymol_script_for_charged_pairs(pdb_file, charged_pairs, output_folder):
    script_filename = os.path.join(output_folder, f"{os.path.basename(pdb_file).replace('.pdb', '')}_highlight.pml")
    
    with open(script_filename, 'w') as script:
        script.write(f"load {pdb_file}\n")
        script.write("show cartoon\n")
        script.write("color green\n")

        for resi1, resn1, resi2, resn2 in charged_pairs:
            script.write(f"select pair_{resi1}_{resi2}, resi {resi1} or resi {resi2}\n")
            script.write(f"color blue, resi {resi1}\n")
            script.write(f"color red, resi {resi2}\n")
            script.write(f"show sticks, pair_{resi1}_{resi2}\n")
        
        script.write("zoom\n")
    print(f"[INFO] PyMOL script saved to {script_filename}")

# Function to process PDB files
def process_pdb_files(pdb_folder, cleaned_pdb_folder, buried_charged_pairs_folder):
    for pdb_file in os.listdir(pdb_folder):
        if pdb_file.endswith(".pdb"):
            pdb_file_path = os.path.join(pdb_folder, pdb_file)
            cleaned_pdb_file = os.path.join(cleaned_pdb_folder, pdb_file)
            
            print(f"Processing {pdb_file}...")
            clean_pdb(pdb_file_path, cleaned_pdb_file)  # Save cleaned PDB to cleaned_pdb_folder
            charged_pairs = calculate_sasa_for_charged_pairs(cleaned_pdb_file)  # Analyze cleaned PDB
            
            if charged_pairs:
                copy_pdb_with_charged_pairs(pdb_file_path, charged_pairs, buried_charged_pairs_folder)  # Save buried charged PDB to separate folder
                generate_pymol_script_for_charged_pairs(cleaned_pdb_file, charged_pairs, buried_charged_pairs_folder)  # Save PyMOL script

# Run the process on all PDB files
process_pdb_files(pdb_folder, cleaned_pdb_folder, buried_charged_pairs_folder)

# Quit PyMOL session
cmd.quit()
