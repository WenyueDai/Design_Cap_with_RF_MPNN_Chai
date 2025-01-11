"""This is a collection of functions useful when dealing with PyMol. Developed
with Python 3.8.5"""

from os.path import isfile
from string import ascii_lowercase, ascii_letters
from glob import glob
from statistics import mean
from random import randint, choice
from re import search

from pymol import cmd, finish_launching
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np

TRIMERS = {"1a": ["01", "02", "05"], "1b": ["15", "16", "17"],
           "2a": ["18", "19", "20"], "2b":["24", "25", "26"],
           "3a": ["27", "28", "29"], "3b": ["21", "22", "23"]}
LINKED_TRIMERS = {"1a": ["02"], "1b": ["05"],
                  "2a": ["08"], "2b": ["14"],
                  "3a": ["15"], "3b": ["11"]}

# A dictionary for translating the three-letter PyMol residue names to the single-letter code.
ONE_LETTER = {"VAL":"V", "ILE":"I", "LEU":"L", "GLU":"E", "GLN":"Q", "ASP":"D", "ASN":"N",
              "HIS":"H", "TRP":"W", "PHE":"F", "TYR":"Y", "GLY":"G", "PRO":"P", "CYS":"C",
              "ARG":"R", "LYS":"K", "SER":"S", "THR":"T", "MET":"M", "ALA":"A"}

####################################################################################################

def add_identifiers():
    """This utility function gives each protein a unique segment identifier and
    potentially a unique chain identifier. Note that a pymol file must already
    be loaded for this function to work.
    """
    cmd.alter("polymer.protein", "segi='00'")
    counter = 0
    while cmd.count_atoms("segi 00") > 0:
        counter += 1
        cmd.alter("bymol first segi 00", f"segi='{(str(counter).zfill(2))}'")
    cmd.sort()
    cmd.dss()

    if counter > 62:
        print("Chain counter is greater than 62, not adding chain identifiers.")
    else:
        chain_names = iter(ascii_letters + "".join(str(list(range(10)))))
        for i in range(counter):
            cmd.alter(f"segi {(str(i + 1).zfill(2))}", f"chain='{next(chain_names)}'")


def mimic_colouring(template_file, template_object, recipient_file, recipient_object):
    """This utility function copies the atom colouring from one pymol .pse
    session to another pymol-compatible file.

    Args:
        template_file (str): path to pymol .pse file with desired colouring
        template_object (str): name of the object in the .pse session that has
            the desired colouring
        recipient_file (str): path to the pymol-compatible file to which the
            colouring should be applied
        recipient_object (str): name of the object in the pymol-compatible file
            to which the desired colouring should be applied

    Returns:
        bool: a statement whether the function has been ran successfully
    """
    for file in [template_file, recipient_file]:
        if not isfile(file):
            print(f"Couldn't find file called {file}, stopping.")
            return False

    finish_launching(["pymol", "-q"])
    cmd.delete("all")
    cmd.load(template_file)
    cmd.remove("hydrogens (not polymer.protein) (not alt ''+A)")
    atoms = pymol_iterate(template_object, dict, "[(model, segi, chain, resi, name)] = color")

    cmd.delete("all")
    cmd.load(recipient_file)
    for atom_identi, color in tqdm(atoms.items()):
        atom_sele = "/".join([recipient_object] + list(atom_identi[1:]))
        if cmd.count_atoms(atom_sele) == 0:
            print(f"Atom {atom_sele} could not be found, continuing.")
        else:
            cmd.color(color, atom_sele)
    return True


def movie_commands():
    """This function creates the frames necessary for a movie showing partial
    (three super-pentamer) capsid expansion for T1 to T3 geometry. Note that a
    suitable pymol file must already be loaded for this function to work.
    """
    cmd.set_view([
        0.103366800,    0.324615300,    0.940180063,
        -0.846459627,    0.525080621,   -0.088233821,
        -0.522317469,   -0.786706686,    0.329053402,
        0.000000000,    0.000000000, -1222.187133789,
        117.797714233,    0.798587799,   43.285903931,
        789.208007812, 1655.166259766,  -20.000000000])
    cmd.color("green", "superpent_one_copy")
    p10_pseud = (cmd.get_coords("p10_pseud") * 1.1).tolist()[0]
    cmd.translate(p10_pseud, "superpent_one_copy", camera=0)
    cmd.rotate(p10_pseud, 36, "superpent_one_copy", camera=0, origin=[0,0,0])
    cmd.color("red", "superpent_two_copy")
    p_10_pseud = (cmd.get_coords("p-10_pseud") * 1.1).tolist()[0]
    cmd.translate(p_10_pseud, "superpent_two_copy", camera=0)
    cmd.rotate(p_10_pseud, 36, "superpent_two_copy", camera=0, origin=[0,0,0])
    cmd.color("blue", "superpent_three_copy")
    _10p_pseud = (cmd.get_coords("10p_pseud") * 1.1).tolist()[0]
    cmd.translate(_10p_pseud, "superpent_three_copy", camera=0)
    cmd.rotate(_10p_pseud, 36, "superpent_three_copy", camera=0, origin=[0,0,0])


def md_radius_change(md_files):
    """This function measures the distance between adjacent trimers (across the
    hexameric interface) over the course of molecular simulation and plots the
    results as a line chart. It can deal both with standard trimers (three
    chains) and linked trimers (one chain).

    Args:
        md_files (list): strings with paths to relevant traj.pdb files
    """
    distances = {}
    for md_file in md_files:
        cmd.delete("all")
        cmd.load(md_file)
        cmd.load_traj(md_file[:-4] + ".dcd", start=2)
        cmd.remove("hydrogens or not polymer.protein")
        add_identifiers()

        dist_over_time = []
        for state in range(1, cmd.count_states() + 1):
            coords = {}
            # 30 proteins in a hexamer with standard trimers, 15 with linked trimers.
            trimers = TRIMERS if len(cmd.get_chains()) == 30 else LINKED_TRIMERS
            for trimer, segi in trimers.items():
                coords[trimer] = cmd.centerofmass(f"name CA & segi {'+'.join(segi)}", state=state)
            mean_dist = mean([euclidean(coords[f"{i}a"], coords[f"{i}b"]) for i in range(1, 4)])
            dist_over_time.append(round(mean_dist, 3))
        distances[int(search(r"_z(\d+)_", md_file).group(1))] = dist_over_time

    for radius in sorted(distances.keys(), reverse=True):
        plt.plot(list(range(1, len(distances[radius]) + 1)), distances[radius])
    plt.legend(sorted(distances.keys(), reverse=True))
    plt.show()


def pymol_iterate(selection, inner_space_type, expression, auxiliary={}): #pylint: disable=dangerous-default-value
    """This utility function decreases the amount of boilerplate code required
    to use the cmd.iterate function.

    Args:
        selection (str): pymol atomic selection to iterate over
        inner_space_type (container): python container (list, dictionary, etc.)
            that will get populated with entries during iteration. Pass it
            without call brackets, i.e. use list not list().
        expression (str): python expression through which entries are added to
            the inner_space container during iteration
        auxiliary (dict, optional): dictionary of additional python objects
            required for the above-described expression to function. Defaults
            to an empty dictionary.

    Returns:
        container: inner_space container filled with data about iterated atoms
    """
    space = {"inner_space": inner_space_type(), **auxiliary}
    cmd.iterate(selection, f"inner_space{expression}", space=space)
    return space["inner_space"]


def renumber(chains, pymol_object="*", resi_offset=0, unique_resi=True):
    """A utility function for renumbering the residues in a chain according to
    their connectivity (i.e., from N terminus to C terminus).

    Args:
        chains (list): which chains' residues to renumber and what order to go
            through the chains in
        pymol_object (str, optional): which object's residues to renumber.
            Defaults to *, which assumes there is only one pymol object in the
            file.
        resi_offset (int, optional): How much to shift the residue numbering.
            Defaults to 0, meaning the renumbered residues will start from 1.
        unique_resi (bool, optional): If multiple chains are specified, the
            numbering for each chain will either continue from where the
            previous chain stopped (if unique_resi is True), or drop back to
            the resi_offset value plus 1 (if unique_resi is False). Defaults to
            True.

    Returns:
        bool: a statement whether the function has been ran successfully
    """
    for chain in chains:
        # Finding the N terminal residue by searching for a nitrogen not bonded to carbonyl carbon.
        msc = f"{pymol_object}/*/{chain}"
        residues = [int(resi) for resi in pymol_iterate(f"{msc}/*/CA", list, ".append(resi)")]
        for resi in residues:
            if cmd.count_atoms(f"{msc}/*/C & bound_to {msc}/{resi}/N") == 0:
                renumbering = [resi]
                residues.remove(resi)
                break
        else:
            print("Couldn't find N-terminal residue, stopping.")
            return False

        # Iteratively finding the next residue. #pylint: disable=undefined-loop-variable
        while len(residues) > 0:
            if cmd.count_atoms(f"{msc}/{resi + 1}/N & bound_to {msc}/{resi}/C") == 1:
                resi += 1
            else:
                for resi_2nd in residues:
                    if cmd.count_atoms(f"{msc}/{resi_2nd}/N & bound_to {msc}/{resi}/C") == 1:
                        resi = resi_2nd
                        break
                else:
                    print(f"Couldn't find the next residue after resi {resi}, stopping.")
                    return False
            renumbering.append(resi)
            residues.remove(resi)

        # Renumbering the residues according to the new scheme.
        for index, resi in enumerate(renumbering):
            cmd.label(f"{msc}/{resi}/*", str(index + 1 + int(resi_offset)))
        cmd.alter(f"{msc}/*/*", "resi=int(label)")
        cmd.sort()
        print(f"Chain {chain} renumbered, continuing.")
        if unique_resi:
            resi_offset += index + 1
    return True


def shorten_sele(resi_list):
    """A utility function for converting a very large list of individual
    residues into groups of consecutive residues. Its main purpose is to
    decrease the number of characters needed to select all these residues in
    PyMol, as that software has a ~3800 character limit per command.

    Args:
        resi_list (list): residue numbers of residues of interest. Both strings
            and integers are accepted.

    Yields:
        str: the first and last residue for a group of consecutive residues,
            separated by a dash. This is the PyMol format for selecting a range
            of residues.
    """
    # Removing duplicates and converting all residues into integers.
    resi_list = sorted(set(int(resi) for resi in resi_list))
    group = resi_list[0:1]
    for resi in resi_list[1:]:
        if resi == group[-1] + 1: # i.e. if the two residues are consecutive
            group.append(resi)
        else:
            yield f"{group[0]}-{group[-1]}" if len(group) > 1 else str(group[0])
            group = [resi]
    yield f"{group[0]}-{group[-1]}" if len(group) > 1 else str(group[0])


def find_pdb_codes(file_type="cif", max_limit=float("inf")):
    """This function generates four-letter strings that could be codes for
    structures in the pdb database and attempts to fetch them. Codes which do
    correspond to a structure are recorded.

    Args:
        file_type (str, optional): What file type to be searching for. Defaults
            to "cif".
        max_limit (int, optional): How many correct pdb codes to find. Defaults
            to float("inf"), meaning positive infinity.

    Yields:
        str: names of correct pdb codes
    """
    found_pdbs, all_characters = [], list(ascii_lowercase) + list(range(10))
    while len(found_pdbs) < max_limit:
        pdb_code = [randint(1, 9)] + [choice(all_characters) for _ in range(3)] # A list
        pdb_code = "".join(str(ch) for ch in pdb_code) # Converting it into one string
        if pdb_code not in found_pdbs:
            cmd.delete("all")
            cmd.fetch(pdb_code, type=file_type, file="temporary_pdb." + file_type)
            if cmd.count_atoms("(all)") > 0:
                found_pdbs.append(pdb_code)
                yield pdb_code


def fetch_sequences(seq_of_interest, search_space, msc=("*", "*", "*"), autoload=True):
    """This utility function searches through a list of pymol-compatible files
    for those that contain a sequence with particular residue numbers and
    residue names.

    Args:
        seq_of_interest (dict): residue numbers (str or int) as keys and
            residue names following the three-leter code (str) as values.
        search_space (list): ABSOLUTE paths to the PyMol-compatile files that
            might contain the sequence of interest
        msc (tuple, optional): the model, segment, and chain names that the
            function should be searching through. Defaults to ("*", "*", "*"),
            meaning everything.
        autoload (bool, optional): whether to display files containing the
            sequence of interest in a PyMol window, one at a time. The search
            gets paused until the user closes the file. Defaults to False.

    Returns:
        dict or False: paths to found files (str) as keys and details of what
            the found sequence was (dict) as values
    """
    for file in search_space:
        if not isfile(file):
            print(f"Couldn't find file called {file}, stopping.")
            return False
    # For unclear reasons, trying to load a PyMol window mid-function doesn't work, so I do it here.
    if autoload:
        finish_launching(["pymol", "-q"])
    msc = f"model {msc[0]} & segi {msc[1]} & chain {msc[2]}"
    # Ensuring correct formatting.
    seq_of_interest = {str(resi): resn for (resi, resn) in seq_of_interest.items()}
    for resi, resn in seq_of_interest.items():
        if resn not in ONE_LETTER:
            print(f"Incorrect residue name {resn} for resi {resi}, stopping.")
            return False

    matching_files = {}
    # Tqdm in combination with autoload prints too much to the terminal.
    for file in (search_space if autoload else tqdm(search_space)): #pylint: disable=superfluous-parens
        cmd.delete("all")
        cmd.load(file)
        selection = f"{msc} & name CA & resi {'+'.join(seq_of_interest.keys())}"
        seq = pymol_iterate(selection, dict, "[resi] = resn")
        if seq == seq_of_interest:
            matching_files[file] = seq
            if autoload:
                cmd.zoom(f"byres ({selection})")
                cmd.show("sticks", "(all)")
                input(f"Press any button to close the file '{file}' and continue searching.")

    if not matching_files:
        print(f"No files containing the sequence '{seq_of_interest}' could be found, stopping.")
    return matching_files


def apply_mutations(mutations, models=r".*", chains="*"):
    """This utility function uses the PyMol protein mutagenesis wizard to
    change residues.

    Args:
        mutations (dict):  residue numbers (str or int) as keys and residue
            names following the three-leter code (str) as values
        models (string, optional): RegEx search pattern for choosing which
            models to mutate. Defaults to r".*", meaning all models.
        chains (string, optional): PyMol selection pattern for choosing which
            chains to mutate. Defaults to "*", meaning all chains.

    Returns:
        bool: a statement whether the function has been ran successfully
    """
    mutations = {str(resi): resn for (resi, resn) in mutations.items()}
    for resi, resn in mutations.items():
        if resn not in ONE_LETTER:
            print(f"Incorrect residue name {resn} for resi {resi}, stopping.")
            return False

    cmd.wizard("mutagenesis")
    #cmd.refresh_wizard() #* Might be needed if multiple wizards are activated.
    wizard = cmd.get_wizard()

    for model in cmd.get_object_list("all"):
        if search(models, model) is None:
            continue
        for chain in cmd.get_chains(f"model {model} & chain {chains}"):
            for resi, resn in mutations.items():
                wizard.set_mode(resn)
                selection = f"{model}//{chain}/{resi}/"
                if pymol_iterate(selection, set, "add(resn)") == resn:
                    print(f"Residue {selection} is already a {resn}, continuing.")
                    continue
                wizard.do_select(selection)
                wizard.apply()
    return True


def calculate_dihedral(points, chiral_center=False):
    """This utility function calculates the dihedral angle between four atoms.
    The atoms are either bonded to each other in a line (A-B-C-D) or form a
    chiral center (with A, C, and D all connecting to B). Maths taken from
    https://stackoverflow.com/questions/20305272/.

    Args:
        points (array): Numpy array containing four rows and three columns,
            corresponding to the x/y/z cartesian coordinates of four atoms
        chiral_center (bool): whether the atoms are bonded in a linear (e.g.
            peptide bond) or branching fashion (e.g. chiral center). Defaults
            to False, meaning linear bonding.

    Returns:
        float or False: dihedral angle value, possible range is from -180 to
            +180.
    """
    try:
        points = np.array(points)
    except: #pylint: disable=bare-except
        print("Points couldn't be converted into an array, stopping.")
        return False
    if points.shape != (4, 3):
        print(f"Incorrect shape of points array {points.shape}, require (4, 3), stopping.")
        return False

    vector_0 = (points[1] - points[0]) * -1
    vector_1 = points[2] - points[1]
    vector_1 /= np.linalg.norm(vector_1) # Normalised vector
    vector_2 = (points[3] - points[2]) if not chiral_center else (points[3] - points[1])

    # Projections of vector_0 or vector_2 onto a plane perpendicular to vector_1.
    projection_0 = vector_0 - np.dot(vector_0, vector_1) * vector_1
    projection_2 = vector_2 - np.dot(vector_2, vector_1) * vector_1

    # Angle between the two projections in a plane is the dihedral angle.
    x_component = np.dot(projection_0, projection_2)
    y_component = np.dot(np.cross(vector_1, projection_0), projection_2)
    return np.degrees(np.arctan2(y_component, x_component))

####################################################################################################

if __name__ == "__main__":
    #pylint: disable=line-too-long
    LOCAL_DIR = "/home/bstka/Work/Coding_playground/"
    #mimic_colouring(LOCAL_DIR + "interface_zones.pse", "contact_map", LOCAL_DIR + "gradual_squish.pse", "100_16")
    #md_radius_change(glob(LOCAL_DIR + "Testing_folder_3/Amber_results/*/traj.pdb"))
    #finish_launching(["pymol", "-q"])
    #apply_mutations({764: "ASN", 1223: "VAL", 1234: "ALA", 1260: "PHE", 761: "GLU", 1468: "GLU",
    #                 1470: "SER", 1720: "GLU", 1737: "LYS", 1723: "ALA"}, r"3kic_tri[A-E]_linked", "a")
    #renumber("abcdefghijklmnopqrstuvwxyz", unique_resi=False)

else: # Loading the script via PyMol
    for func in (add_identifiers, movie_commands, pymol_iterate, renumber, shorten_sele,
                 apply_mutations, calculate_dihedral):
        cmd.extend(func.__name__, func)
