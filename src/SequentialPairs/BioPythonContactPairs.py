#!/Users/ameliaantrim/anaconda/bin/python

'''
Parses a PDB file. 

Input:
(1) File: PDB file containing the structure to be parsed
(2) Cutoff: The maximum distance (in angstroms) between two residues for them to be 
    considered "in contact"

Class Attributes:
(1) self.file = the PDB file containing the protein information
(2) self.structure = the parsed PDB file containing the Bio.PDB.structure object
(3) self.Pairs = contact pairs from the contact map and structure info
(4) self.NamePairs
(5) self.ContactPairs = contact pairs on different chains

'''

import sys
import math
import os
from Bio.PDB import *
from collections import defaultdict

class PDBParse:
    def __init__(self, File, cutoff):
        self.file = File
        print 'HI WE ARE USING THE CORRECT PARSER'

        # Initialize the Bio.PDB parser and get the structure from PDB file
        parser = PDBParser(PERMISSIVE=1)
        self.structure = parser.get_structure("struct", File)
        
        # Initialize the arrays containing the pair information
        self.ContactPairs = []
        self.NamePairs = []    
        self.Pairs = []
        
        # Dictionary to convert amino acid 3-letter names to 1-letter names for
        # conversion to other applications and string formats
        AAConverter = { "LEU": "L",
        "GLY": "G", "LYS": "K", "SER": "S", "VAL": "V", "ARG": "R", "THR": "T", 
        "PRO": "P", "ILE": "I", "MET": "M","PHE": "F","TYR": "Y","CYS": "C",
        "CYX": "C","CYD": "C", "TRP": "W","HID": "H","HIE": "H","HIS": "H","HIP":"H",
        "ALA": "A", "GLU": "E", "GLN": "Q", "ASP": "D", "ASN": "N"}

        #get contact pairs
        #iterate through the pairs of residue numbers
        # Compare each chain
        for residue1 in self.structure.get_residues():
            for residue2 in self.structure.get_residues():
                # If no carbon alpha, can't get the coordinates; likely water
                if 'CA' not in residue1:
                    continue
                elif 'CA' not in residue2:
                    continue
                # If the distance between them is less than the cutoff distance
                elif (residue1['CA'] - residue2['CA']) <= cutoff:
                
                    # Convert to the one-letter names if not already converted
                    if len(residue1.get_resname()) > 1:
                        residue1.resname = AAConverter[residue1.get_resname()]
                    if len(residue2.get_resname()) > 1:
                        residue2.resname = AAConverter[residue2.get_resname()]
                    # Add to pairs and name pairs
                    self.Pairs.append((residue1, residue2))
                    self.NamePairs.append((residue1.get_resname(), residue2.get_resname()))
                    if residue1.get_parent() != residue1.get_parent():
                        self.ContactPairs.append((residue1, residue2))
