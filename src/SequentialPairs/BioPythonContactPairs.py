#!/Users/ameliaantrim/anaconda/bin/python

import sys
import math
import os
from Bio.PDB import *
from collections import defaultdict

'''
Parses a PDB file. 
Contents:
	self.file = the PDB file containing the protein information
	self.Length = the number of amino acids in the chain
	self.seq = amino acid sequence of chain
	self.StartingResidue = first residue in chain
	self.ContactMap = contact residues, defined as the distance between
		them being no more than "cutoff" angstroms (here, 8)
	self.DistanceMap = the distances between contact pair residues
	self.Pairs = contact pairs from the contact map and structure info
	self.Distance = distances between the above contact pairs

Warning: Only works on a single chain of the protein. A discontinuous protein warning
often suggests that the protein in question has multiple chains. 
'''
class PDBParse:
    def __init__(self, File):
    	self.file = File
        self.Coord = []
        ResNumList =[]
        parser = PDBParser(PERMISSIVE=1)
        self.structure = parser.get_structure("struct", File)
        self.Pairs = []
        #AAConverter = defaultdict(char)
        AAConverter = { "LEU": "L",
        "GLY": "G",
        "LYS": "K",
        "SER": "S",
        "VAL": "V",
        "ARG": "R",
        "THR": "T",
        "PRO": "P",
        "ILE": "I",
        "MET": "M","PHE": "F","TYR": "Y","CYS": "C","CYX": "C","CYD": "C",\
        "TRP": "W","HID": "H","HIE": "H","HIS": "H","HIP":"H"}
        AAConverter["ALA"] = "A"
        AAConverter["GLU"] = "E"
        AAConverter["GLN"] = "Q"
        AAConverter["ASP"] = "D"
        AAConverter["ASN"] = "N"
        self.seq = []
        # NOTE: B IS FOR BAD RESIDUE NAMES (NAG)
        for residue in self.structure.get_residues():
        	# Add the single-letter version of the amino acid name
        	self.seq.append(residue.get_resname())
        	print residue.get_resname()
        	
        #self.Length = len(self.seq)
        #self.StartingResidue = ResNumList[0]

		#get contact pairs
        cutoff = 8
        self.ContactMap = []
        self.DistanceMap = []
        self.NamePairs = []
        #iterate through the pairs of residue numbers
        # Compare each chain
        
        # CONTACT PAIRS
        for chain1 in self.structure.get_chains():
        	print 'CHAIN' + chain1
			# To all the other chains
			for chain2 in self.structure.get_chains():
				# Compare all residues in each
				for residue1 in chain1.get_residues():
					# To all residues in all other chains
					for residue2 in chain2.get_residues():
						# If no carbon alpha, can't get the coordinates
						if 'CA' not in residue1:
							continue
						elif 'CA' not in residue2:
							continue
						# If they are in the same chain, ignore
						elif residue1.get_parent() == residue2.get_parent():
							if residue1.get_resname() in AAConverter.keys() and residue2.get_resname in AAConverter.keys():
								self.Pairs.append((residue1, residue2))
								self.NamePairs.append((AAConverter[residue1], AAConverter[residue2]))
							#print (residue1, residue2)
							continue
						# Check the distance between them and check if it's below the cutoff (here, 8)
						elif (residue1['CA'] - residue2['CA']) <= cutoff:
							# If they are in contact, append them to the list of interface contact pairs
							self.ContactPairs.append((residue1, residue2))
							# The sphere dict tells how many times that sphere appears in a pair that is a contact pair
    	        			if residue1.get_resname() in AAConverter.keys() and residue2.get_resname in AAConverter.keys():
	    	        			self.contact_name_pairs.append((AAConverter[residue1.get_resname()], AAConverter[residue2.get_resname()]))
    	        				self.Pairs.append((residue1, residue2))
    	        				self.NamePairs.append((AAConverter[residue1], AAConverter[residue2]))
    	        		else:
    	        			continue
        
        '''for residue1 in self.structure.get_residues():
        	for residue2 in self.structure.get_residues():
        		# Control for NAG residues
        		if not "CA" in residue1:
        			continue
        		if not "CA" in residue2:
        			continue
        		ca1 = residue1['CA']
        		ca2 = residue2['CA']
	            #self.ContactMap.append(CM_line)
    	        #self.DistanceMap.append(DM_line)
    	        #check the distance
    	        if ca1-ca2 <= cutoff:
    	        	# If they are in contact, append them to the list of interface contact pairs
    	        	self.Pairs.append((residue1, residue2))
    	        	print self.Pairs[0]
    	        	if residue1.get_resname() in AAConverter.keys() and residue2.get_resname in AAConverter.keys():
	    	        	self.name_pairs.append((AAConverter[residue1.get_resname()], AAConverter[residue2.get_resname()]))'''
    	        			

		

'''Parse and represent a single line in a PDB file'''
class LineParse:
    def __init__(self, line):
        self.AtomName = line[12:16].strip()
        self.ResName = line[17:20].strip()
        self.ResNum = int(line[22:26])
        self.X = float(line[30:38])
        self.Y = float(line[38:46])
        self.Z = float(line[46:54])