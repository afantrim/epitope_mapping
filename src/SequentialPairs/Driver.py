#!/Users/ameliaantrim/anaconda/bin/python

'''
Driver class in which I can enter files and parameters.

Pipeline:
(1) Checks validity of inputs
(2) Calls SequentialPairs to turn the peptides into a dictionary of pairs with the
    values being the number of times that pair appears in the peptide dataset.
    **check in on valued spheres***
(3) Calls PyMolColor to generate two PyMol batch files to color both the lines and 
    the spheres
    ***check spheres***
(4) Generates real contact pairs if that information is available
    ***maybe think about making this a separate script

To do: 
- Make a list of the real positives
- Make a list of the predicted positives
- Add to the script functionality to calculate the true and false positives
- Fix the reduce function
- Decide between PDBParse/Biopython contact pairs and remove the other 
- Ditto w/ Biopython_sequential_pairs and PDBParse (probably delete PDBParse and then
    rename Biopython_sequntial_pairs to PDBParse because this is a more logical name

'''

from collections import defaultdict
import sys
import math
import os
from ASAParse import ASAParse
from BioPythonContactPairs import PDBParse
from SequentialPairs import SequentialPairs
from SubGroups import SubGroups
from PyMolColor import PyMolColor
from Bio.PDB import *
from PyMolColorSpheres import ColorResidues
from Statistics import Statistics
from PrecisionRecall import PrecisionRecall


'''Check whether a file exists and can be opened'''
def check_file(the_file):
    if not os.path.isfile(the_file):
        print the_file + " does not exist."
        return False
    elif not os.access(the_file, os.R_OK):
        print "Problem opening file " + the_file
        return False
    return True

'''
Driver to run this utility from the command line
Argument format: [protein file from PDB] [surface residue file from ASA] [sequence file]

'''
def main():
    if len(sys.argv) >= 4: # Get files from command line
        peptide_file, asa_file, sequence_file = sys.argv[1], sys.argv[2], sys.argv[3]
        if len(sys.argv) >= 5:
            if sys.argv[4] == 0:
                reduced = False # Whether we will used the reduced alphabet
            else:
                reduced = True
        else: # Default is to not be reduced
            reduced = False
    else:
        opt = raw_input("Wrong # of file inputs. Enter names manually (y/n)? ")
        if opt[0] == "y":
            peptide_file = raw_input("Enter path to peptide file: ")
            asa_file = raw_input("Enter path to ASA file: ")
            sequence_file = raw_input("Enter path to sequence file: ")
        else:
            print "Exiting. Re-try with new input (see README).\n"
            exit(0)
    cutoff = 17
    
    # Create the antibody/peptide model and color the pairs in pymol
    if not check_file(asa_file) and check_file(peptide_file) and check_file(sequence_file):
        print "poop"
        exit()

    # Parse out the sequential pairs in the peptide dataset
    model = SequentialPairs(peptide_file, reduced)

    # Get a dictionary representation of the surface area
    asa_dict = ASAParse(asa_file).SA_dict

    # Parse the PDB file to get antigen structural and sequence information 
    pdb_info = PDBParse(sequence_file, cutoff)

    # Make a PyMOL batch file to color the residues according to whether they 
    # are estimated to be part of the epitope
    spheres = ColorResidues(model.pair_dict, sequence_file, pdb_info, asa_dict, reduced)
    spheres.create_residues_file()

    # Get information for co-crystallized structure
    #antibody_asa_file = raw_input("Enter the antibody ASA file: ")
    #antibody_pdb_file = raw_input("Enter the antibody PDB file: ")
    #co_pdb_file = raw_input("Enter co-crystallized PDB file: ")
    antibody_asa_file = sys.argv[5]
    antibody_pdb_file = sys.argv[6]
    co_pdb_file = sys.argv[7]

    # Model the antibody and the antigen in co-crystallized form 
    contact_pairs = PDBParse(co_pdb_file, cutoff).ContactPairs

    # Then make a CSV file to compare the false and true positives
    positives = Statistics(spheres.positives, contact_pairs)

    plot = PrecisionRecall()
    plot.add_threshhold(positives.true_positives, positives.false_positives)
    
if __name__ == "__main__":
    main()

