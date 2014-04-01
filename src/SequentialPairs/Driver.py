#!/Users/ameliaantrim/anaconda/bin/python

from collections import defaultdict
import sys
import math
import os
from ASAParse import ASAParse
from PDBParse import PDBParse
from SequentialPairs import SequentialPairs
from SubGroups import SubGroups
from PyMolColor import PyMolColor
from ContactPairs import ContactPairs
from Bio.PDB import *

'''Outputs a PyMOL batch file to color amino acids as spheres based on their occurence
in the peptide dataset '''
class ResidueColor:
	def __init__(self):
		pass
	# Get data from SequentialPairs
	# For each pair, take the sum of the lines between that and its contact pair
	# NOTE: Work on this as soon as you are finished converting to Biopython. Make
	# things cleaner and automate 

'''Check whether a file exists and can be opened'''
def check_file(the_file):
	if not os.path.isfile(the_file):
		print the_file + " does not exist."
		return False
	elif not os.access(the_file, os.R_OK):
		print "Problem opening file " + the_file
		return False
	return True

'''Driver to run this utility from the command line
Argument format: [protein file from PDB] [surface residue file from ASA] [sequence file]
'''
def main():
	if len(sys.argv) >= 4: # Get files from command line
		peptide_file, asa_file, sequence_file = sys.argv[1], sys.argv[2], sys.argv[3]
		if len(sys.argv) >=5:
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
	
	# Create the antibody/peptide model and color the pairs in pymol
	if check_file(asa_file) and check_file(peptide_file) and check_file(sequence_file):
		model = SequentialPairs(peptide_file, asa_file, sequence_file, reduced)
		#print (model.pair_dict.keys())
		#print model.pair_dict.values()
		colors = PyMolColor(model.pair_dict, model.valued_spheres, sequence_file, asa_file, reduced)
		colors.create_colors_file()
	else:
		print "Problem with files. Exiting.\n"
		exit()
	
	# MAKE THIS OPTIONAL
	# antigen_ASA_file, antibody_ASA_file, antigen_pdb_file, antibody_pdb_file, co_PDB_file	
	antibody_asa_file = raw_input("Enter the antibody ASA file: ")
	antibody_pdb_file = raw_input("Enter the antibody PDB file: ")
	co_pdb_file = raw_input("Enter co-crystallized PDB file: ")
	contact_pairs = ContactPairs(asa_file, antibody_asa_file, sequence_file, antibody_pdb_file, co_pdb_file)
	

if __name__ == "__main__":
	main()













