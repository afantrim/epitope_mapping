#!/Users/ameliaantrim/anaconda/bin/python


import os
import sys
sys.path.append('/Users/ameliaantrim/anaconda')
from ASAParse import ASAParse
from PDBParse import PDBParse
from Bio.PDB import *

'''
Takes two molecules (usually an antibody and an antigen) in complex and the ASA
files of each molecule giving surface exposure values. From this information,
it provides a list of the residue types and numbers of the contact pairs between
the two molecules.
'''
class ContactPairs:
	def __init__(self, antigen_ASA_file, antibody_ASA_file, antigen_pdb_file, antibody_pdb_file, co_PDB_file):
		# Get the PDB parser from BioPython
		parser = PDBParser()
		
		# Parse the separate antibody and antigen PDB files
		self.antigen_pdb_info = parser.get_structure("antigen", antigen_pdb_file)
		self.antibody_pdb_info = parser.get_structure("antibody", antibody_pdb_file)
		
		# Get combined PDB files with co-crystallized structure
		self.crystal_structure = parser.get_structure("co-crystallized", co_PDB_file)
		self.crystal_pairs = PDBParse(co_PDB_file).Pairsc
		
		# Parse the ASA surface area files for array of surface areas
		self.antigen_asa = ASAParse(antigen_ASA_file).SA_dict
		self.antibody_asa = ASAParse(antibody_ASA_file).SA_dict
		
		'''for pair in self.co_pdb_info.Pairs:
			print pair'''
		self.get_contact_pairs()
		self.display_contact_pairs()
	
	'''
	For surface residue A in antigen: 
		For surface residue B in antibody:
			If (A, B) or (B, A) is a pair in co_pdb_info:
				Add this pair to the contact pairs
	'''
	def get_contact_pairs(self):
		# Start with the antigen surface b/c it's the smaller dataset; less to deal with
		self.potential_residues = []
		# Get all the necessary pairs
		for residue1 in self.antigen_structure.get_residues():
			for residue2 in self.antibody_structure.get_residues():
				if (residue1, residue2) in self.crystal_pairs:
					self.potential_residues.append((residue1, residue2))
	'''
	Prints out all the contact pairs
	'''
	def display_contact_pairs(self):
		for pair in self.potential_residues:
			print pair
			print "\n"
	
	'''
	Makes a file that will highlight all the contact pairs
	'''
	def make_file(self):
		# Set up file; this will be only for the co-crystallized file
		# Highlight contact pairs in complex
		with open("contact_pairs.txt", "w") as conctact_file:
			for pair in self.potential_residues:
				print "color residue " + str(pair[0])
				print "color residue " + str(pair[1])
					
def main():
	antigen_asa = sys.argv[1]
	antibody_asa = sys.argv[2]
	antigen_pdb = sys.argv[3]
	antibody_pdb = sys.argv[4]
	co_pdb = sys.argv[5]
	ContactPairs(antigen_asa, antibody_asa, antigen_pdb, antibody_pdb, co_pdb)

if __name__ == "__main__":
	main()