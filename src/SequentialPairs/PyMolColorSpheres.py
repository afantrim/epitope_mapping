#!/Users/ameliaantrim/anaconda/bin/python

import math
import os
from collections import defaultdict
from ASAParse import ASAParse
from BioPythonContactPairs import PDBParse
from Biopython_Sequential_Pairs import SequentialPairs
from SubGroups import SubGroups
from Bio.PDB import *

MIN_ST = 30

'''Create PyMOL files to color the molecule according to number of surface pairs'''
class ColorResidues:
        #get the pair_dict from SequentialPairs
    def __init__(self, pair_dict, pdb_file, pdb_info, asa_dict, reduced):
        # Make info accessible to all functions
        self.pdb_info = pdb_info
        self.pair_dict = pair_dict
        self.pdb_file = pdb_file
        self.colors = []
        self.colors.append("myred, [1.00 , 0.00 , 0.00]")
        self.colors.append("myorange, [1.00, 0.50, 0.00]")
        self.colors.append("myyellow, [0.95, 0.78, 0.00]")
        self.colors.append("myblue, [0.02, 0.50, 0.72]")
        self.colors.append("myvioletpurple, [0.55, 0.25, 0.60]")

        # Translate each residue to the reduced amino acid alphabet if applicable
        if reduced:
            for residue in self.pdb_info.structure.get_residues():
                SubGroups.residue_reduce(residue)

        # Parse the ASA file to get the surface accessibility d
        self.asa_dict = asa_dict
        self.count_residues()

    def get_residue_color(self, residue_value):
        global MIN_ST
        print residue_value
        value = float(residue_value - MIN_ST)
        maximum = float(max(self.residues.itervalues())-MIN_ST)
        x = (value / maximum) * 4
        print value
        print maximum
        print x
        return self.colors[int(x)].split()[0]

    def count_residues(self):
        print "hi"
        self.residues = defaultdict(int)
        for pair in self.pdb_info.Pairs:
            if (pair[0].get_resname(), pair[1].get_resname()) in self.pair_dict.keys():
                self.residues[pair[0].id[1]] += self.pair_dict[(pair[0].get_resname(), pair[1].get_resname())]
                self.residues[pair[1].id[1]] += self.pair_dict[(pair[0].get_resname(), pair[1].get_resname())]

    def create_residues_file(self):
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

            #dictionary to hold all the pairs to be colored
        self.max_residue_color = self.get_max_residue()
        self.residue_dict = defaultdict(int)
        self.positives = []

        for residue in self.residues.keys():
            if self.residues[residue] < MIN_ST:
                continue
            pymol_file.write("sel ABC = (resi " + str(residue) + " and n. CA)\n")
            pymol_file.write("show spheres, ABC\n")
            pymol_file.write("color " + self.get_residue_color(self.residues[residue]) + "ABC\n")
            pymol_file.write("deselect\n\n")
            # Keep a running tally of the residues estimated to be in these pairs
            self.positives.append(residue)
        
        pymol_file.close()
        log_file.close()

    def get_max_residue(self):
        return max(self.residues.itervalues())
