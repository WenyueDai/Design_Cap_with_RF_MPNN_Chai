import os
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis import Universe
from MDAnalysis.analysis import rms
from scipy.stats import gaussian_kde

"""
conda activate getcontact
Use after RFdiffusion.
Filter cap based on criteria:
- Distance: how far cap is away from main body
- Angle: angle between cap and main body
- RMSD: how similar between cap and main body
"""

# Set default font size for all plots
plt.rcParams.update({'font.size': 32})

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
    sequence_string = "".join(sequence.values())
    return sequence_string

# Function to determine glycine-based N cap and C cap lengths (RFdiffusion generate all G)
def determine_glycine_caps(sequence):
    n_cap_len, c_cap_len = 0, 0
    for i, res in enumerate(sequence):
        if res == 'G':
            n_cap_len += 1
        else:
            break
    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] == 'G':
            c_cap_len += 1
        else:
            break
    return n_cap_len, c_cap_len

# Distance and angle calculations
def calculate_angle(point_a, point_b, point_c):
    vec_a = point_a - point_b
    vec_c = point_c - point_b
    cosine_angle = np.dot(vec_a, vec_c) / (np.linalg.norm(vec_a) * np.linalg.norm(vec_c))
    angle = np.degrees(np.arccos(np.clip(cosine_angle, -1.0, 1.0)))
    return angle

# Function to calculate RMSD, distance, and angle for N and C caps
def calculate_cap_parameters(pdb_file):
    sequence = parse_pdb_sequence(pdb_file, chain_id='A')
    if not sequence:
        print(f"No valid sequence in {pdb_file}. Skipping.")
        return None

    n_cap_len, c_cap_len = determine_glycine_caps(sequence)
    if n_cap_len == 0 or c_cap_len == 0:
        print(f"Skipping {pdb_file} due to insufficient N or C cap length.")
        return None
    
    u = Universe(pdb_file)

    # N cap and C cap atoms
    n_cap = u.select_atoms(f"resid 1:{n_cap_len}")
    c_cap = u.select_atoms(f"resid {len(sequence) - c_cap_len + 1}:{len(sequence)}")
    main_body_n = u.select_atoms(f"resid {n_cap_len + 1}:{n_cap_len * 2}")
    main_body_c = u.select_atoms(f"resid {len(sequence) - 2 * c_cap_len + 1}:{len(sequence) - c_cap_len}")

    if len(n_cap) != len(main_body_n) or len(c_cap) != len(main_body_c):
        return None

    # RMSD calculations
    rmsd_n_analysis = rms.RMSD(n_cap, main_body_n).run()
    rmsd_c_analysis = rms.RMSD(c_cap, main_body_c).run()
    rmsd_n = rmsd_n_analysis.results.rmsd[-1, 2]
    rmsd_c = rmsd_c_analysis.results.rmsd[-1, 2]

    n_cap_com = n_cap.center_of_mass()
    c_cap_com = c_cap.center_of_mass()
    center_line_com = (main_body_n.center_of_mass() + main_body_c.center_of_mass()) / 2.0
    n_cap_distance = np.linalg.norm(n_cap_com - center_line_com)
    c_cap_distance = np.linalg.norm(c_cap_com - center_line_com)
    n_cap_angle = calculate_angle(n_cap_com, center_line_com, main_body_n.center_of_mass())
    c_cap_angle = calculate_angle(c_cap_com, center_line_com, main_body_c.center_of_mass())

    return {
        'pdb_file': pdb_file,
        'n_cap_rmsd': rmsd_n,
        'n_cap_distance': n_cap_distance,
        'n_cap_angle': n_cap_angle,
        'c_cap_rmsd': rmsd_c,
        'c_cap_distance': c_cap_distance,
        'c_cap_angle': c_cap_angle
    }

# Function to save parameters to a text file
def save_parameters_to_file(results, output_file):
    with open(output_file, 'w') as f:
        # Write header
        f.write("PDB_File\tN_Cap_RMSD\tN_Cap_Distance\tN_Cap_Angle\tC_Cap_RMSD\tC_Cap_Distance\tC_Cap_Angle\n")
        # Write each result line by line
        for res in results:
            f.write(f"{res['pdb_file']}\t{res['n_cap_rmsd']:.3f}\t{res['n_cap_distance']:.3f}\t{res['n_cap_angle']:.2f}\t"
                    f"{res['c_cap_rmsd']:.3f}\t{res['c_cap_distance']:.3f}\t{res['c_cap_angle']:.2f}\n")

# Function to plot individual histograms for each parameter with preset bins
def plot_individual_histograms(data, title_prefix, output_folder):
    parameters = ['RMSD', 'Distance', 'Angle']
    axis_limits = {
        'RMSD': (0, 20),
        'Distance': (0, 30),
        'Angle': (0, 60)
    }
    
    for i, param in enumerate(parameters):
        values = data[:, i]
        
        # Use fixed bin settings for consistency across histograms
        n_bins = bin_settings[param]
        
        plt.figure(figsize=(10, 8))
        plt.hist(values, bins=n_bins, color='orange')
        plt.title(f"{title_prefix} {param} Distribution", pad=20)
        plt.xlabel(param)
        plt.ylabel("Frequency")
        
        # Set axis limits
        plt.xlim(axis_limits[param][0], axis_limits[param][1])
        plt.ylim(0, 100)  # Cap frequency to ensure readability
        
        filename = f"{title_prefix.lower()}_{param.lower()}_distribution.png"
        plt.savefig(os.path.join(output_folder, filename), bbox_inches='tight')
        plt.close()

# Function to plot pairwise parameter comparisons with density coloring and fixed axis limits
def plot_pairwise_parameters(data, title_prefix, output_folder):
    pairs = [
        ('RMSD', 'Distance', data[:, 0], data[:, 1], (0, 20, 0, 30)),
        ('RMSD', 'Angle', data[:, 0], data[:, 2], (0, 20, 0, 60)),
        ('Distance', 'Angle', data[:, 1], data[:, 2], (0, 30, 0, 60))
    ]
    
    for x_label, y_label, x, y, limits in pairs:
        # Calculate density for coloring
        xy = np.vstack([x, y])
        density = gaussian_kde(xy)(xy)
        
        plt.figure(figsize=(10, 8))
        scatter = plt.scatter(x, y, c=density, cmap='viridis', marker='o')
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(f"{title_prefix} {x_label} vs {y_label}")
        
        # Set axis limits
        plt.xlim(limits[0], limits[1])
        plt.ylim(limits[2], limits[3])
        
        # Add color bar for density
        cbar = plt.colorbar(scatter)
        cbar.set_label('Density')
        
        # Save plot
        filename = f"{title_prefix.lower()}_{x_label.lower()}_vs_{y_label.lower()}.png"
        plt.savefig(os.path.join(output_folder, filename), bbox_inches='tight')
        plt.close()

# Main function for Script 1
def calculate_caps_parameters(input_folder, output_folder, parameter_file):
    results = []
    for pdb_file in os.listdir(input_folder):
        if pdb_file.endswith('.pdb'):
            pdb_path = os.path.join(input_folder, pdb_file)
            parameters = calculate_cap_parameters(pdb_path)
            if parameters:
                results.append(parameters)
    
    # Convert results to numpy arrays for easier plotting
    n_cap_data = np.array([[r['n_cap_rmsd'], r['n_cap_distance'], r['n_cap_angle']] for r in results])
    c_cap_data = np.array([[r['c_cap_rmsd'], r['c_cap_distance'], r['c_cap_angle']] for r in results])
    
    # Generate individual histograms for N Cap and C Cap data
    plot_individual_histograms(n_cap_data, "N Cap", output_folder)
    plot_individual_histograms(c_cap_data, "C Cap", output_folder)
    
    # Generate pairwise plots for N Cap and C Cap data
    plot_pairwise_parameters(n_cap_data, "N Cap", output_folder)
    plot_pairwise_parameters(c_cap_data, "C Cap", output_folder)
    
    # Save parameters to file for use in the second script
    save_parameters_to_file(results, parameter_file)

# Example usage
# Predefined bin counts based on parameter range
bin_settings = {
    'RMSD': 18,      # Adjust based on data range (0–20) and manual (increase if frequency>100, decrease if too many bin)
    'Distance': 15,  # Adjust based on data range (0–30) and manual
    'Angle': 18      # Adjust based on data range (0–60) and manual
}
input_folder = "/home/eva/0_bury_charged_pair/5_Pipeline_r1_extenlen_25_40/step1_rfdiffusion_m8/output/pdb"
output_folder = "/home/eva/0_bury_charged_pair/5_Pipeline_r1_extenlen_25_40/step1_rfdiffusion_m8/1_cap_calculation"
os.makedirs(output_folder, exist_ok=True)
parameter_file = os.path.join(output_folder, "cap_parameters.txt")
calculate_caps_parameters(input_folder, output_folder, parameter_file)
