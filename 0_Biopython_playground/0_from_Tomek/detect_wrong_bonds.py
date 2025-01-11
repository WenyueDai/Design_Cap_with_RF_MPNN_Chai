"""This script detects and possibly fixes incorrect bonds in pymol-compatible
files. It can be run both on Linux and Windows. Developed with Python 3.8.5"""

from collections import defaultdict
from os.path import isfile, dirname
from os import chdir
import pickle

from pymol import cmd, finish_launching

from pymol_utilities import pymol_iterate, calculate_dihedral

####################################################################################################

def create_bond_rules(input_file):
    """This function detects all bonds in a pymol-compatible file and compiles
    this information into a rule-set of permitted bonds, which it saves as a
    .pickle file. Likewise, you must be certain that the input file contains no
    incorrect bonds, and has all 20 natural amino-acids present. Note that you
    won't be able to open the .pickle file in a text editor to view it - use
    the view_ruleset function instead.

    Args:
        input_file (str): path to the pymol-compatible input file

    Returns:
        output_file or False: string specifying the path to the output .pickle
            file
    """
    if not isfile(input_file):
        print(f"Couldn't find file called {input_file}, stopping.")
        return False
    cmd.delete("all")
    cmd.load(input_file)
    cmd.remove("hydrogens (not polymer.protein) (not alt ''+A)")

    # The data structure here is: sets inside dictionaries inside a dictionary
    allowed_bonds = defaultdict(lambda: defaultdict(set))
    atoms = pymol_iterate("(all)", dict, "[(ID, model)] = [name, resn, resi]")
    for atom_key, (name, resn, resi) in atoms.items():

        cmd.select(f"bound_to ({atom_key[1]} & id {atom_key[0]})")
        neighbor_set = set()
        neighbors = pymol_iterate(r"%sele", list, ".append((name, resi))")
        for (neighbor_name, neighbor_resi) in neighbors:
            # neighbor_resi is re-purposed here to express whether the two bonded atoms are in the
            # same residue (0) or adjacent residues (-1 or 1). An exception is made for disulphide
            # sulphurs, where the two residues can be any distance from each other, sequence-wise.
            if neighbor_name == name == "SG":
                neighbor_resi = "any"
            else:
                neighbor_resi = int(neighbor_resi) - int(resi)
            neighbor_set.add((neighbor_name, neighbor_resi))

        allowed_bonds[resn][name].union(neighbor_set)

    # Manually adding the special OXT oxygen atom (found on the C-terminus) to each residue.
    for residue_name in allowed_bonds:
        allowed_bonds[residue_name]["C"].add(("OXT", 0))
        # Making a one-member set containing a two-member tuple.
        allowed_bonds[residue_name]["OXT"] = set([("C", 0)])

    # Changing the file type to .pickle but keeping the rest of the filename intact.
    output_file = ".".join(input_file.split(".")[:-1]) + ".pickle"
    with open(output_file, "wb") as file:
        pickle.dump(allowed_bonds, file) # Converting the allowed_bonds dictionary to a file.
    return output_file


def view_ruleset(ruleset_file):
    """This utility function allows you to view and potentially edit the
    .pickle file specifying permitted bonds.

    Args:
        ruleset_file (str): path to the .pickle file with bond rules

    Returns:
        bool: a statement whether the function has been ran successfully
    """
    if not isfile(ruleset_file):
        print(f"Couldn't find file called {ruleset_file}, stopping.")
        return False
    with open(ruleset_file, "r+b") as file:
        ruleset_dict = pickle.load(file)
        #* Carry out any edits over here.
        for residue_name, residue_atoms in ruleset_dict.items():
            print("RESIDUE ", residue_name) # residue_name is a string
            for atom, neighbor_atoms in residue_atoms.items():
                # atom is a string, neighbor_atoms is a set of two-member tuples
                print("  ", atom, neighbor_atoms)
        file.seek(0)
        pickle.dump(ruleset_dict, file)
        file.truncate()
    return True


def error_detection(pymol_file, ruleset_file, fix=False, objects_of_interest=None):
    """This function goes through all the atoms in a pymol-compatible file and
    compares their bonds to a user-specified bond ruleset. If the ruleset was
    well-constructed, all possible wrong bonds will be detected and (if you set
    fix=True) automatically removed.

    Args:
        pymol_file (str): path to the pymol-compatible input file
        ruleset_file (str): path to the .pickle file with bond rules
        fix (bool): whether the function should remove wrong bonds that it
            finds. Defaults to False.
        objects_of_interest (list, optional): list of strings, specifying which
            pymol objects the function should analyse for wrong bonds. Defaults
            to None, which means that all objects in the file are analysed.

    Returns:
        set: set of two-member tuples, which specify the atom ID and object
            name for each atom with wrong bonds.
    """
    for file_path in [pymol_file, ruleset_file]:
        if not isfile(file_path):
            print(f"Couldn't find file called {file_path}, stopping.")
            return False
    cmd.delete("all")
    cmd.load(pymol_file)
    objects_of_interest = "(all)" if objects_of_interest is None else " ".join(objects_of_interest)
    # Working with hydrogens is currently too complicated for this script.
    cmd.select(f"polymer.protein & ({objects_of_interest}) & not hydrogens")
    with open(ruleset_file, "rb") as file:
        ruleset_dict = pickle.load(file)

    wrong_bonds, key_errors = set(), set()
    atoms = pymol_iterate(r"%sele", dict, "[(ID, model)] = [name, resn, resi, chain]")
    for atom_key, (name, resn, resi, chain) in atoms.items():

        cmd.select(f"polymer.protein & not hydrogens & bound_to ({atom_key[1]} & id {atom_key[0]})")
        for neighbor in pymol_iterate(r"%sele", list, ".append([name, resi, chain, ID, model])"):

            neighbor[1] = "any" if neighbor[0] == name == "SG" else int(neighbor[1]) - int(resi)
            try:
                # This is the main check for wrong bonds.
                if (neighbor[0], neighbor[1]) not in ruleset_dict[resn][name]:
                    wrong_bonds.update([atom_key, (neighbor[3], neighbor[4])])
                    if fix:
                        cmd.unbond(f"{atom_key[1]} & id {atom_key[0]}",
                                   f"{neighbor[4]} & id {neighbor[3]}")
            except KeyError:
                key_errors.add((resn, name))
            # This check looks for the rare case of an incorrect bond being made between two chains
            # by residues with the same or +1 / -1 residue number (which might not be caught above).
            if neighbor[2] != chain and not neighbor[0] == name == "SG":
                wrong_bonds.update([atom_key, (neighbor[3], neighbor[4])])

    print("The wrong bonds were:", wrong_bonds)
    print("The key errors were:", key_errors) # This should be empty.
    return wrong_bonds


def detect_cis_backbone(pymol_file, states=(1, 2, 1)):
    """This auxiliary function measure the dihedral angles of all peptide bonds
    in a pymol-compatible file and categorises them as cis/trans/neither.

    Args:
        pymol_file (str): path to the pymol-compatible input file
        states (tuple, optional): which states to analyse. Uses the same
            arguments as the built-in range function. Defaults to (1, 2, 1),
            meaning only state 1 is considered.

    Returns:
        bool: a statement whether the function has been ran successfully
    """
    if not isfile(pymol_file):
        print(f"Couldn't find file called {pymol_file}, stopping.")
        return False
    cmd.delete("all")
    cmd.load(pymol_file)
    cmd.color("white", "(all)")
    cmd.show_as("spheres", "polymer.protein & name N")
    cmd.color("green", "rep spheres")

    not_trans = defaultdict(set)
    backbone = "polymer.protein & name N+CA+C"
    identifiers = pymol_iterate(backbone, list, ".append((model, segi, chain, resi, resn))")

    for state in range(*states):
        coords = cmd.get_coords(backbone, state=state) # Same order as with pymol_iterate
        if coords is None:
            print(f"State {state} does not exist, continuing.")
            continue
        zipped = list(zip(coords, identifiers))
        # Grouping the atoms into peptide bonds (CA, C, N, CA).
        peptide_bonds = [zipped[i:i+4] for i in range(1, len(zipped)-4, 3)]

        for atoms in peptide_bonds:
            points, identi = zip(*atoms)
            # All the atoms need to be from the same model, segi, and chain to form a peptide bond.
            if len(set(item[:3] for item in identi)) != 1:
                continue
            angle = abs(round(calculate_dihedral(points)))
            for cutoff, category in [(130, "trans"), (50, "neither"), ( 0, "cis")]:
                if angle >= cutoff:
                    if category != "trans":
                        not_trans[category].add(identi[2])
                        print(f"{identi[2]} state {state} was {category}, with an angle of {angle}")
                    break

    for category in sorted(not_trans, reverse=True): # neither, then cis
        for identi in not_trans[category]:
            cmd.color("yellow" if category == "cis" else "red", "/".join(identi[:-1]) + "/N")
    return True


def detect_aminoacid_chirality(pymol_file, states=(1, 2, 1)):
    """This auxiliary function measure the chirality at the alpha carbon of all
    non-glycine aminoacids in a pymol-compatible file and categorises them as L
    (by far the predominant form in nature) or D.

    Args:
        pymol_file (str): path to the pymol-compatible input file
        states (tuple, optional): which states to analyse. Uses the same
            arguments as the built-in range function. Defaults to (1, 2, 1),
            meaning only state 1 is considered.

    Returns:
        bool: a statement whether the function has been ran successfully
    """
    if not isfile(pymol_file):
        print(f"Couldn't find file called {pymol_file}, stopping.")
        return False
    cmd.delete("all")
    cmd.load(pymol_file)
    cmd.color("white", "(all)")
    cmd.show_as("spheres", "polymer.protein & name CA & not resn GLY")
    cmd.color("green", "rep spheres")

    l_aminoacids = set()
    backbone = "polymer.protein & name N+CA+C+CB & not resn GLY"
    identifiers = pymol_iterate(backbone, list, ".append((model, segi, chain, resi, resn))")

    for state in range(*states):
        coords = cmd.get_coords(backbone, state=state) # Same order as with pymol_iterate
        if coords is None:
            print(f"State {state} does not exist, continuing.")
            continue
        zipped = list(zip(coords, identifiers))
        # Grouping the atoms into chiral groups (N, CA, C, CB).
        chiral_centers = [zipped[i:i+4] for i in range(0, len(zipped), 4)]

        for atoms in chiral_centers:
            points, identi = zip(*atoms)
            angle = round(calculate_dihedral(points, chiral_center=True))
            if angle > 0:
                l_aminoacids.add(identi[2])
                print(f"{identi[2]} state {state} had D chirality, with an angle of {angle}")

    for identi in l_aminoacids:
        cmd.color("red", "/".join(identi[:-1]) + "/CA")
    return True

####################################################################################################

if __name__ == "__main__":
    chdir(dirname(__file__))
    finish_launching(["pymol", "-qx"])
    detect_cis_backbone("D:\\Schoolwork\\Freeline_postdoc\\AAV_codebase\\Testing_folder_6\\2022-11-14_4918_10mut_T3_hexamer_hairpin1_sim\\cis.pdb", (1, 74, 2))
    #detect_aminoacid_chirality("D:\\Schoolwork\\Freeline_postdoc\\AAV_codebase\\Testing_folder_5\\3kic_triC_linked_T3_10mut_z103_capsid_amber_sim\\no_cis_T3.pdb", (1, 110, 5))
    exit()
    #pickled_ruleset = create_bond_rules("6fwt.pdb")
    view_ruleset("6fwt.pickle")
    cmd.delete("all")
    incorrect_bonds = error_detection("6ih9_triE_T3_morphed.cif", "6fwt.pickle", False, ["tri_0"])
    cmd.select("wrong_bonds", " ".join([f"{item[1]} & id {item[0]}" for item in incorrect_bonds]))
    cmd.show("spheres", "wrong_bonds")
