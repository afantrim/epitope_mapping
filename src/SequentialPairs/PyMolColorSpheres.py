#!/Users/student/anaconda/bin/python

'''
.. module:: PyMolColorSpheres
   :platform: Unix, Windows
   :synopsis: Outputs batch files to color residues according to different paramters. 

.. moduleauthor:: Amelia F. Antrim <amelia.f.antrim@gmail.com>

'''

import math
import sys
import os
from collections import defaultdict
from ASAParse import ASAParse
from PDBParse import PDBParse
from Biopython_Sequential_Pairs import SequentialPairs
from SubGroups import SubGroups
from Bio.PDB import *
import PDBUtil
from MatrixReduce import *
from Blosum import *

MIN_ST = -100

class ColorResidues:

    '''
    .. function: __init__(self, pair_dict, pdb_file, pdb_info, asa_dict, reduced=False, scoring_matrix=None)
    
    Create PyMOL files to color the molecule according to number of surface pairs.
    
    :param pair_dict: Residue number pairs mapped to number of times they appear in peptide dataset.
    :type pair_dict: dictionary (key=int, value=int or float)
    :param pdb_file: PDB file containing the structure of the protein to be colored with batch file.
    :type pdb_file: str(filename)
    :param pdb_info: Parsed structure of protein to be colored with batch file.
    :type pdb_info: `PDBParse <file:///Users/student/Dropbox/epitope_mapping/src/SequentialPairs/docs/_build/html/index.html#module-PDBParse>`_ instance.
    :param asa_dict: dictionary of residue numbers mapped to their relative accessible surface area
    :type asa_dict: dictionary (key=int, value=float)
    :param reduced: boolean indicating whether we are using a reduced alphabet
    :type reduced: bool
    :param scoring_matrix: matrix mapping pairs of amino acids to a substitutability index
    :type scoring_matrix: nested dictionaries (key=int,value=dictionary(key=int,value=float or int))
    
    '''
    def __init__(self, pair_dict, pdb_file, pdb_info, asa_dict, reduced=False, scoring_matrix=Blosum.CSERZO):
        # Make info accessible to all functions
        self.pdb_info = pdb_info
        self.pair_dict = pair_dict
        self.pdb_file = pdb_file
        self.asa_dict = asa_dict
        self.colors = []
        self.colors.append("myred, [1.00 , 0.00 , 0.00]")
        self.colors.append("myorange, [1.00, 0.50, 0.00]")
        self.colors.append("myyellow, [0.95, 0.78, 0.00]")
        self.colors.append("myblue, [0.02, 0.50, 0.72]")
        self.colors.append("myvioletpurple, [0.55, 0.25, 0.60]")

        # Translate each residue to the reduced amino acid alphabet if applicable
        '''if reduced:
            for residue in self.pdb_info.structure.get_residues():
                SubGroups.residue_reduce(residue)
            self.count_reduced_residues()

        else:
            self.count_residues(scoring_matrix)'''
        self.count_residues(scoring_matrix)

    def _get_residue_color(self, residue_value):
        '''Calculates how this residue should be colored based on its frequency in the peptide 
        sequential pairs, relative to other residue pairs.
        '''
        global MIN_ST
        value = float(residue_value - MIN_ST)
        maximum = float(max(self.residues.itervalues())-MIN_ST)
        x = (value / maximum) * 4
        return self.colors[int(x)].split()[0]

    def count_residues(self, scoring_matrix):
        '''
        Iterates through pairs and determines how many times each contact pair appears in 
        the peptide dataset; for when using a scoring matrix.
        
        '''
        self.residues = defaultdict(float)
        pairs = []
        for pair in self.pdb_info.Pairs:
            print pair
            pairs.append(pair)
        while len(pairs) > 0:
            pair = pairs.pop()
            for true_pair in self.pair_dict.keys():
                score = MatrixReduce.score_pair(true_pair, pair, scoring_matrix)
                if score > 0:
                    self.residues[pair[0].id[1]] += score*(self.pair_dict[true_pair]**2)
                    #print score*self.pair_dict[true_pair]
                    self.residues[pair[1].id[1]] += score*(self.pair_dict[true_pair]**2)
                    #print pair, true_pair
        print self.residues

    def count_reduced_residues(self):
        '''
        Iterates through pairs and determines how many times each contact pair appears in 
        the peptide dataset; this is for when using a reduced amino acid alphabet rather than the
        scoring matrix.
        
        '''
        for pair in self.pdb_info.Pairs:
            if (pair[0].get_resname(), pair[1].get_resname()) in self.pair_dict.keys():
                if pair[0].id[1] not in self.residues.keys():
                    self.residues[pair[0].id[1]] = 1
                else:
                    self.residues[pair[0].id[1]] += self.pair_dict[(pair[0].get_resname(), pair[1].get_resname())]
                if pair[1].id[1] not in self.residues.keys():
                    self.residues[pair[1].id[1]] = 1
                else:
                    self.residues[pair[1].id[1]] += self.pair_dict[(pair[0].get_resname(), pair[1].get_resname())]

    def create_residues_file(self):
        '''
    	Generates a PDB batch file to color residues of a protein based on the frequency of each
    	residue's appearance in a pair represented in the peptide dataset, and how many times that 
    	pair appears in the peptide dataset.

    	'''
    
        global MIN_ST

        print "Creating residue coloring file..."
        # Create a pymol batch file for commands
        pymol_file = open(self.pdb_file[:-4]+"_residues.txt", "w")
        log_file = open(self.pdb_file[:-4]+"_residues.csv", "w")

        # Start writing the header info
        pymol_file.write("hide lines\n")
        pymol_file.write("set dash_width, 2\n")
        pymol_file.write("set dash_gap, 0\n\n")

        # Write in colors
        for color in self.colors:
            pymol_file.write("set_color " + color + "\n")
        pymol_file.write("\n")

        # Dictionary to hold all the pairs to be colored
        self.max_residue_color = self.get_max_residue()
        self.residue_dict = defaultdict(int)
        self.positives = []
        
        # Iterate through pairs to create file
        for residue in self.residues.keys():
            # Must be above the minimum occurrence threshhold
            if self.residues[residue] < MIN_ST:
                continue
            else:
                pymol_file.write("sel ABC = (resi " + str(residue) + " and n. CA)\n")
                pymol_file.write("show spheres, ABC\n")
                pymol_file.write("color " + self._get_residue_color(self.residues[residue]) + "ABC\n")
                pymol_file.write("deselect\n\n")
                # Keep a running tally of the residues estimated to be in these pairs
                self.positives.append(residue)
                #print self.positives
        
        pymol_file.close()
        log_file.close()

    def get_max_residue(self):
        '''
        Determines the maximum index of a residue's appearances in pairs in the peptide dataset.
        
        :returns: int -- maximum number of times a residue appears in a pair in peptides
        '''
        return max(self.residues.itervalues())
