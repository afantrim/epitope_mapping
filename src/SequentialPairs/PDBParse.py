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


'''
class PDBParse:
    def __init__(self, File, cutoff):
    	self.file = File
        self.Coord = []
        ResNumList =[]
        parser = PDBParser(PERMISSIVE=1)
        self.structure = parser.get_structure("struct", File)
        self.Pairs = []
        #AAConverter = defaultdict(char)
        AAConverter = { "LEU": "L", "GLY": "G", "LYS": "K", "SER": "S",
        "VAL": "V", "ARG": "R", "THR": "T", "PRO": "P", "ILE": "I",
        "MET": "M","PHE": "F","TYR": "Y","CYS": "C","CYX": "C","CYD": "C",\
        "TRP": "W","HID": "H","HIE": "H","HIS": "H","HIP":"H", 
        "ASN" : "N", "ALA" : "A", "GLU": "E", "GLN" : "Q", "ASP" : "D"}
        
        # seq will contain the names of all the residues in single-letter format
        self.seq = []
        for residue in self.structure.get_residues():
        	self.seq += AAConverter[residue.get_resname()]
        
        # The length and first residue of the sequence
        self.Length = len(self.seq)
        self.StartingResidue = self.seq[0]

		#get contact pairs
        self.ContactMap = []
        self.DistanceMap = []
        self.NamePairs = []
        #iterate through the pairs of residue numbers
        # Compare each chain
        
        # Compare each residue
        for residue1 in self.structure.get_residues():
        	# To each other residue
        	for residue2 in self.structure.get_residues():
        	
        		# Control for water molecules
        		if not is_aa(residue1):
        			continue
        		if not is_aa(residue2):
        			continue
        		
        		# Control for NAG residues
        		atom1 = 'CA'
        		atom2 = 'CA'
        		
        		# If there is no alpha carbon, use the beta carbon to track distance (glycine)
        		if not residue1.has_id('CA'):
        			atom1 = 'CB'
        			
        		if not residue2.has_id('CA'):
        			atom2 = 'CB'
        		
        		# Get the distance between the two residues
        		ca1 = residue1[atom1]
        		ca2 = residue2[atom2]
        		distance = (ca1 - ca2)
	            
    	        # Check the distance
    	        if distance <= cutoff:
    	        	print distance
    	        	# If they are in contact, append them to the list of interface contact pairs
    	        	self.Pairs.append((residue1, residue2))
					# Append the residue names into a separate pair for easier access
    	        	if residue1.get_resname() in AAConverter.keys() and residue2.get_resname() in AAConverter.keys():
	    	        	self.NamePairs.append((AAConverter[residue1.get_resname()], AAConverter[residue2.get_resname()]))
    	        			
