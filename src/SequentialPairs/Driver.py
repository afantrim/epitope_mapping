#!/Users/student/anaconda/bin/python

'''
.. module:: Driver
   :platform: Unix, Windows
   :synopsis: A sample driver class to analyze an epitope. 

.. moduleauthor:: Amelia F. Antrim <amelia.f.antrim@gmail.com>

'''

from collections import defaultdict
import sys
import math
import os
from ASAParse import ASAParse
from PDBParse import PDBParse
from SequentialPairs import SequentialPairs
from SubGroups import SubGroups
#from PyMolColor import PyMolColor
from Bio.PDB import *
from PyMolColorSpheres import ColorResidues
from Statistics import Statistics
from PrecisionRecall import PrecisionRecall
from PDBUtil import PDBUtil
import numpy
from TestThreshhold import TestThreshhold


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
            
    cutoff = 12.0
    
    # Create the antibody/peptide model and color the pairs in pymol
    if not check_file(asa_file) and check_file(peptide_file) and check_file(sequence_file):
        print "poop"
        exit()

    # Parse out the sequential pairs in the peptide dataset
    model = SequentialPairs(peptide_file, reduced)

    # Get a dictionary representation of the surface area
    antigen_asa_dict = ASAParse(asa_file).SA_dict
    
    # Parse the PDB file to get antigen structural and sequence information 
    pdb_info = PDBParse(sequence_file, cutoff, antigen_asa_dict)
    pdb_info.make_pairs()

    # Make a PyMOL batch file to color the residues according to whether they 
    # are estimated to be part of the epitope
    spheres = ColorResidues(model.pair_dict, sequence_file, pdb_info, antigen_asa_dict, reduced)
    spheres.create_residues_file()

    # Get information for co-crystallized structure
    antibody_asa_file = sys.argv[5]
    antibody_pdb_file = sys.argv[6]
    co_pdb_file = sys.argv[7]
    
    # Get the accessible surface area of the antibody
    antibody_asa_dict = ASAParse(antibody_asa_file).SA_dict
    
    # Model the antibody and the antigen in co-crystallized form 
    contact_structure = PDBParse(co_pdb_file, asa_dict=antibody_asa_dict)
    
    # The first structure in the model (there is only one, but the organization of the PDB requires this)
    contact_model = contact_structure.structure

    # Make the plots
    TestThreshhold.make_plot(spheres.residues, contact_model, contact_model["C"], dist=12.0)
    
    # Then make a CSV file to compare the false and true positives
    #positives = Statistics(spheres.positives, contact_pairs)
    
if __name__ == "__main__":
    main()

