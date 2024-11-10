import os
import subprocess
import re
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.exceptions import SelectionError
import pandas as pd
import matplotlib.pyplot as plt

# Set global font size for plots
plt.rcParams.update({'font.size': 16})

# Define the paths
chai_folder = "/home/eva/0_bury_charged_pair/0_With_initial_pdb/3_Chai/1lxa"
rdfiffusion_pdb = "/home/eva/0_bury_charged_pair/0_With_initial_pdb/1_interesting_start/1lxa_HDHD/1lxa_HDHD.pdb"  # Single PDB file to compare
analysis_folder = "/home/eva/0_bury_charged_pair/0_With_initial_pdb/4_Analysis"
rmsd_output_file = os.path.join(analysis_folder, "rmsd_results.txt")
scatter_plot_file = os.path.join(analysis_folder, "rmsd_vs_chai_score.png")
rmsd_vs_chai_output_file = os.path.join(analysis_folder, "rmsd_vs_chai_score.txt")
summary_file = os.path.join(chai_folder, "best_prediction_summary.txt")

# Lists to store RMSD and Chai score results for plotting and sorting
rmsd_data = []
chai_scores = []
sorted_data = []

# Function to calculate global RMSD for Cα atoms only
def calculate_global_rmsd(structure1, structure2):
    try:
        u1 = mda.Universe(structure1)
        u2 = mda.Universe(structure2)
        alignment = rms.RMSD(u1, u2, select="name CA")
        alignment.run()
        return alignment.results.rmsd[-1, 2]
    except SelectionError as e:
        print(f"SelectionError: {e}")
        return None

# Function to convert .cif to .pdb using PyMOL
def convert_cif_to_pdb(cif_file, output_folder):
    pdb_file = os.path.join(output_folder, os.path.splitext(os.path.basename(cif_file))[0] + ".pdb")
    try:
        # Call PyMOL for conversion
        pymol_cmd = f"pymol -c -d 'load {cif_file}; save {pdb_file}; quit'"
        subprocess.run(pymol_cmd, shell=True, check=True)
        return pdb_file
    except subprocess.CalledProcessError as e:
        print(f"Error converting {cif_file} to PDB using PyMOL: {e}")
        return None

# Function to read Chai scores from the summary file
def read_chai_scores(summary_file):
    chai_scores = {}
    with open(summary_file, "r") as f:
        lines = re.findall(r"(sequence_\d+\.cif):\s+([\d\.]+)", f.read())
        for file_name, score in lines:
            chai_scores[file_name] = float(score)
    return chai_scores

# Main processing loop
chai_scores_dict = read_chai_scores(summary_file)

with open(rmsd_output_file, "w") as rmsd_file, open(rmsd_vs_chai_output_file, "w") as rmsd_vs_chai_file:
    rmsd_file.write("File Name\tGlobal RMSD (Å)\tChai Score\n")
    rmsd_vs_chai_file.write("File Name\tGlobal RMSD (Å)\tChai Score\n")

    for chai_file in os.listdir(chai_folder):
        if chai_file.endswith(".cif"):
            chai_file_path = os.path.join(chai_folder, chai_file)

            # Convert .cif to .pdb for chai_file
            converted_pdb_path = convert_cif_to_pdb(chai_file_path, chai_folder)
            if converted_pdb_path:
                global_rmsd = calculate_global_rmsd(rdfiffusion_pdb, converted_pdb_path)
                chai_score = chai_scores_dict.get(chai_file, None)

                if global_rmsd is not None and chai_score is not None:
                    rmsd_data.append(global_rmsd)
                    chai_scores.append(chai_score)

                    rmsd_file.write(f"{chai_file}\t{global_rmsd:.4f}\t{chai_score:.4f}\n")
                    rmsd_vs_chai_file.write(f"{chai_file}\t{global_rmsd:.4f}\t{chai_score:.4f}\n")
                    sorted_data.append((chai_file, global_rmsd, chai_score))

# Create scatter plot of RMSD vs. Chai score
if rmsd_data and chai_scores:
    plt.figure(figsize=(10, 6))
    plt.scatter(rmsd_data, chai_scores, alpha=0.7, color='b')
    plt.xlabel("Global RMSD (Å)")
    plt.ylabel("Chai Score")
    plt.title("RMSD vs. Chai Score")
    plt.tight_layout()
    plt.savefig(scatter_plot_file, format="png", dpi=300)
    plt.show()
    print(f"Scatter plot saved as {scatter_plot_file}.")
else:
    print("No data available for scatter plot.")

# Sort and save RMSD and Chai score data
if sorted_data:
    df_sorted = pd.DataFrame(sorted_data, columns=["file_name", "global_rmsd", "chai_score"])
    df_sorted = df_sorted[(df_sorted['global_rmsd'] <= 10) & (df_sorted['chai_score'] >= 0.12)]

    # Find Pareto optimal points
    pareto_front = []
    for index, row in df_sorted.iterrows():
        is_dominated = False
        for _, candidate in df_sorted.iterrows():
            if (candidate['global_rmsd'] < row['global_rmsd'] and candidate['chai_score'] >= row['chai_score']) or \
               (candidate['global_rmsd'] <= row['global_rmsd'] and candidate['chai_score'] > row['chai_score']):
                is_dominated = True
                break
        if not is_dominated:
            pareto_front.append(row)

    pareto_df = pd.DataFrame(pareto_front)
    sorted_output_file = os.path.join(analysis_folder, "sorted_rmsd_chai_scores.txt")
    pareto_df.to_csv(sorted_output_file, index=False, sep='\t')
    print(f"Sorted RMSD and Chai score data saved to {sorted_output_file}.")
