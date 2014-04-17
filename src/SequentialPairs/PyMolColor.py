#!/Users/ameliaantrim/anaconda/bin/python

'''
Creates PyMol files to color molecules according to the number of times each contact
pair appears in the peptide dataset. 

Inputs:
(1) pair_dict: A dictionary mapping each amino acid pair in the peptide dataset to the
    number of times it appears in that dataset
(2) pdb_file: The PDB file containing the structural and sequence information of the 
    antigen
(3) asa_file: The file containing the accessible surface area information
    ***Is this used?***
    ***Think about passing this in as an already-parsed dictionary or removing
    altogether by filtering earlier***
(4) reduced: A boolean indicating whether to use the given reduced amino acid alphabet
    ***May have to change this once multiple alphabets are incorporated***

Outputs:
(1) A file containing information to color the lines indicating pairs
(2) A file containing information to color residues based on the sum of pairs around them

'''

import sys
import math
import os
from ASAParse import ASAParse
from PDBParse import PDBParse
from SequentialPairs import SequentialPairs
from SubGroups import SubGroups
from Bio.PDB import *

MIN_ST = 2

'''Create PyMOL files to color the molecule according to number of surface pairs'''
class PyMolColor:
	#get the pair_dict from SequentialPairs
	def __init__(self, pair_dict, pdb_info, asa_dict, reduced):
		# Make info accessible to all functions
		self.pdb_info = pdb_info
		self.pair_dict = pair_dict
		print pair_dict
		
		# Translate each residue to the reduced amino acid alphabet if applicable
		if reduced:
			for residue in self.pdb_info.structure.get_residues():
				SubGroups.residue_reduce(residue)
		
		# Define the PyMOL colors to be used in the file
		# NOTE: Look into automating this via the "spectrum" function after converting to Biopython
		self.colors = ["myred, [1.00 , 0.00 , 0.00]", "myorange, [1.00, 0.50, 0.00]", \
		"myyellow, [0.95, 0.78, 0.00]", "myblue, [0.02, 0.50, 0.72]",\
		 "myvioletpurple, [0.55, 0.25, 0.60]"]
	
	'''Creates PyMOL batch file for coloring residues and drawing colored links'''
	def create_colors_file(self):
		#create a pymol batch file for commands
		global MIN_ST
		
		
		pymol_file = open(self.pdb_info.file[:-4]+".txt", "w")
		print "Creating line coloring file..."
		
		#start writing the header info
		pymol_file.write("show cartoon\n")
		pymol_file.write("cartoon loop\n")
		pymol_file.write("hide lines\n")
		pymol_file.write("set dash_width, 2\n")
		pymol_file.write("set dash_gap, 0\n\n")
		
		#write in colors
		for color in self.colors:
			pymol_file.write("set_color " + color + "\n")
		pymol_file.write("\n")
			
		#dictionary to hold all the pairs to be colored
		self.max_color = self.get_max()
		
		i = 0
		#get the range of values to be colored
		#set the bin values within that range
		#color the whole thing gray
		for pair in self.pdb_info.Pairs:
			#get the residues at each pair number
			#if represented in peptide data
			#set color according to number of times it appears
			#create a line between the alpha carbons of the two residues
			#NOTE: am getting weird int pairs; fix this when youhave time
			name_pair = (pair[0].get_resname(), pair[1].get_resname())
			if name_pair not in self.pair_dict.keys() or self.pair_dict[name_pair] < MIN_ST:
			    continue
			
			x = str(pair[0].id[1])
			y = str(pair[1].id[1])
			print x
			print y
			print self.pair_dict[name_pair]
			
			this_dist = "mydist_" + str(pair[0].id[1]) + "_" + str(pair[1].id[1])
			pymol_file.write("distance " + this_dist + ", ")
			pymol_file.write("(resi " + x + " and n. CA), ")
			pymol_file.write("(resi " + y + " and n. CA)\n")
			pymol_file.write("color " + self.get_color(name_pair).split()[0] \
					+ " " + this_dist + "\n")

			#use a new line to separate it out
			pymol_file.write("\n")
			i += 1
				
		#write the closing info
		pymol_file.write("hide label\n")
		pymol_file.write("set opaque_background, off\n")
		pymol_file.write("set ray_shadows, off\n")
		pymol_file.write("ray 2048, 1536\n")
		print "Done."
		#pymol_file.write("png saveFilename.png\n")

	'''Puts the color of a single reside into a bin in rainbow array'''
	def get_color(self, pair):
		#take the value of the pair
		#and compare it to (int)(value/max_value)*6
		#then the array is [r][o][y][g][b][v]
		#change -1 to -minimum_color
		pair_value = self.pair_dict[pair]
		x = ((float(pair_value-(MIN_ST-1)) / float(self.max_color-(MIN_ST-1))) * 4)
		return self.colors[int(x)]
		
	'''Returns the maximum number of pairs between any amino acid types in pair_dict'''
	def get_max(self):
		#print self.pair_dict
		return max(self.pair_dict.itervalues())
