# Design capping motif for beta-solenoid protein

This repository provides Python 3 scripts for designing capping motifs for beta-solenoid proteins using a pipeline of RFdiffusion, ProteinMPNN, and Chai-1.

Designed N cap from top view:
<img src="/images/N_cap_top.png" alt="N_cap_top" width=400 height=200>

Designed N cap from side view:
<img src="/images/N_cap_side.png" alt="N_cap_side" width=400 height=200>

Designed C cap from top view:
<img src="/images/C_cap_top.png" alt="C_cap_top" width=400 height=200>

Designed C cap from side view:
<img src="/images//C_cap_side.png" alt="C_cap_side" width=400 height=200>


### Workflow:
1. **Structure Repacking**: The beta-solenoid structure (PDB: 3ult) is repacked with valine at the N- and C-termini and glutamate in the main body.
   
2. **Capping Motif Generation**: RFdiffusion generates capping motifs ranging from 25–40 amino acids, applying block adjacency and secondary structure guide scaffolds.

3. **Initial Filtering**: Capping motifs are evaluated for structural feasibility based on RMSD, distance, and angle between the cap and main body. Protein structures are visualized by randomly selecting 10 candidates.

4. **Separation and Design**: The N- and C-terminal caps are separated, and the main body retains its original structure. RFdiffusion is applied once to design both caps.

5. **Sequence Addition**: PyRosetta restores the wild-type sequence to the main body (PDB: 3ult). SolubleMPNN designs the capping motifs and integrates them with a strand of the wild-type main body.

6. **Structure Prediction**: Chai-1 predicts the structure of the new beta-solenoid. The RMSD between Chai-1 and RFdiffusion structures is compared, and Chai-1 aggregate scores are plotted to determine cutoff thresholds for both RMSD and aggregate scores.

7. **Final Filtering**: Structures are filtered based on:
   - Hydrogen bond formation between the cap and main body.
   - Buried surface area of the cap (calculated as: Main body + Cap - Whole structure).
   - Percentage of helices in the capping motif.
   - Number of atoms in close proximity to the cap.

**Key Criteria**: The percentage of helices and buried surface area are critical for beta-solenoid applications.

## Setup

The script was developed with biopython 1.83, matplotlib 3.8.3, MDAnalysis 2.7.0, numpy 1.26.4, pandas 2.2.0, pymol 2.5.8, pyrosetta 2024.42+release.3366cf78a3, scipy 1.12.0, torch 2.4.1.post100. It has not been tested for version compatibility, so it’s recommended to create a virtual environment with these specific versions to avoid issues. 

## Usage

It is recommend to go through each script named from step0 to step7, and provide the path for the necessary pdb file to run.

## Contributing

Currently not open to contributions, but comments with suggestions are welcome. This repository is being maintained by [Wenyue Dai](https://github.com/WenyueDai).

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/TomaszKaminski-netizen/capsid-expansion/blob/master/LICENSE.txt) file for details.
