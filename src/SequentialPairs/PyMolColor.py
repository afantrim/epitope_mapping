#!/Users/ameliaantrim/anaconda/bin/python

from collections import defaultdict
import sys
import math
import os
from ASAParse import ASAParse
from PDBParse import PDBParse
from SequentialPairs import SequentialPairs
from SubGroups import SubGroups
from Bio.PDB import *

'''Create PyMOL files to color the molecule according to number of surface pairs'''
class PyMolColor:
	#get the pair_dict from SequentialPairs
	def __init__(self, pair_dict, valued_spheres, pdb_file, asa_file, reduced):
		# Make info accessible to all functions
		self.pdb_file = pdb_file
		self.pair_dict = pair_dict
		
		# Get pure PDB structure information 
		parser = PDBParser(PERMISSIVE=1)
		self.pdb_struct = parser.get_structure("self", pdb_file)
		
		# Parse the PDB file via the self-created parser to extract additional information 
		self.pdb_info = PDBParse(pdb_file)
		
		# Translate each residue to the reduced amino acid alphabet if applicable
		self.reduced = reduced
		if reduced:
			self.pdb_info.seq = []
			for residue in self.pdb_info.seq:
				# String reduce the name of the residue (this is the only part that will change)
				self.pdb_info.seq += SubGroups.str_reduce(residue) 
			for residue in self.pdb_struct.get_residues():
				SubGroups.str_reduce(residue.get_resname())
		#os.chdir("/Users/ameliaantrim/Dropbox/epitope_mapping/results")
		
		#for pair in self.pdb_info.Pairs:
		#	print pair[0]
		self.valued_spheres = valued_spheres
		
		# Parse the ASA file to get the surface accessibility data
		self.asa_values = ASAParse(asa_file).SA
		
		# Define the PyMOL colors to be used in the file
		# NOTE: Look into automating this via the "spectrum" function after converting to Biopython
		self.colors = ["myred, [1.00 , 0.00 , 0.00]", "myorange, [1.00, 0.50, 0.00]", \
		"myyellow, [0.95, 0.78, 0.00]", "mygreen, [0, 0.53, 0.22]", "myblue, [0.02, 0.50, 0.72]",\
		 "myvioletpurple, [0.55, 0.25, 0.60]"]
		#print self.pdb_info.Pairs
		
		# Make the file to show spheres
		self.create_spheres_file()
	
	
	'''Creates PyMOL batch file for coloring residues and drawing colored links'''
	def create_colors_file(self):
		#create a pymol batch file for commands
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
		
		#get the range of values to be colored
		#set the bin values within that range
		#color the whole thing gray
		i = 0 
		print "LENGTH OF PAIRS:"
		print len(self.pdb_info.Pairs)
		print "LENGTH OF NAME PAIRS: "
		print len(self.pdb_info.NamePairs)
		for pair in self.pdb_info.Pairs:
			#get the residues at each pair number
			#if represented in peptide data
			#set color according to number of times it appears
			#create a line between the alpha carbons of the two residues
			#NOTE: am getting weird int pairs; fix this when youhave time
			if type(pair[0]) == "int":
				print "weird int residue"
				continue
			this_dist = "mydist_" + str(pair[0].id[1]) + "_" + str(pair[1].id[1])
			pymol_file.write("distance " + this_dist + ", ")
			pymol_file.write("(resi " + str(pair[0].id[1]) + " and n. CA), ")
			pymol_file.write("(resi " + str(pair[1].id[1]) + " and n. CA)\n")
			print "PRE-CALL PAIR"
			print pair
			pymol_file.write("color " + self.get_color(self.pair_dict[self.pdb_info.NamePairs[i]]).split()[0] \
					+ " " + this_dist + "\n")
			#look through lineparse (pdb_info.Coord); get Coord with corresponding ResNum
			#use the pair number to index into the sequence
			
			#use a new line to separate it out
			pymol_file.write("\n")
			i+=1 # Move along in name pairs list
				
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
		x = 0
		if pair in self.pair_dict.keys():
			pair_value = pair_dict[self.name_pair(pair)]
			x = ((float(pair_value-1) / float(self.max_color-1)) * 5)
		return self.colors[int(x)]
		
	'''Returns the maximum number of pairs between any amino acid types in pair_dict'''
	def get_max(self):
		#print self.pair_dict
		return max(self.pair_dict.itervalues())
	
	'''Returns the maximum number of occurrences of pairs containing a single 
	amino acid residue'''
	def get_max_sphere(self):
		return max(self.valued_spheres.itervalues())
	
	'''Gets the number of occurrences of each individual residue included in a pair
	from the peptide file'''
	def create_spheres_file(self):
		#create a pymol batch file for commands
		pymol_file = open(self.pdb_file[:-4]+"_spheres.txt", "w")
		
		#start writing the header info
		pymol_file.write("hide lines\n")
		pymol_file.write("set dash_width, 2\n")
		pymol_file.write("set dash_gap, 0\n\n")
		
		#write in colors
		for color in self.colors:
			pymol_file.write("set_color " + color + "\n")
		pymol_file.write("\n")
			
		#dictionary to hold all the pairs to be colored
		self.max_sphere_color = self.get_max_sphere()
		
		self.sphere_dict = defaultdict(int)
		for pair in self.pdb_info.Pairs:
			# Create a sphere_dict where each residue ID is a key and
			# when you come across it in a pair add the number of times that pair
			# appears in the peptide dataset
			residue1 = pair[0].id[1]
			residue2 = pair[1].id[1]
			if pair in self.pair_dict:
				self.sphere_dict[residue1]+= self.pair_dict[pair]
				self.sphere_dict[residue2] += self.pair_dict[pair]
		
		for sphere in self.sphere_dict.keys():
			pymol_file.write("sel ABC = (resi " + sphere + "and n. CA)\n")
			pymol_file.write("show spheres, ABC\n")
			pymol_file.write("color " + get_sphere_color(sphere_dict[sphere]) + "ABC\n")
			pymol_file.write("deselect\n\n")
			
		pymol_file.close()
				
		'''sel ABC = (resi 100 and n. CA)
		show spheres, ABC
		color red, ABC <- use the color value that corresponds
	 	w/ number of times in peptide dataset
		deselect'''
	

	def get_sphere_color(self, sphere_value):
		print str(sphere_value) + "\n"
		x = ((float(sphere_value-1) / float(self.max_sphere_color-1)) * 5)
		return self.colors[int(x)].split()[0]

			