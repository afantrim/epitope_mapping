#!/Users/ameliaantrim/anaconda/bin/python

import sys
import math
import os
from collections import defaultdict
from ASAParse import ASAParse
from SubGroups import SubGroups

'''Represents a set of sequential pairs of amino acids derived from a set of peptide
binding data, maps them to the surface pairs of the antibody'''
class SequentialPairs:
	def __init__(self, peptide_file, asa_file, sequence_file, reduced):
		self.sequence_file = sequence_file
		self.asa_file = asa_file
		self.peptide_file = peptide_file
		with open(peptide_file, "r") as file:
			self.peptides = file.read().split('\n')
		self.surface_residues = ASAParse(self.asa_file).SA_dict
		#self.reduced_surface = SubGroups(self.surface_residues).reduced_string
		self.valued_spheres = defaultdict(int) # to be filled in in make_pairs()
		if reduced:
			print "SURFACE RESIDUES: \n"
			print self.surface_residues
			temp = self.peptides
			self.peptides = []
			for peptide in temp:
				self.peptides.append(SubGroups.str_reduce(peptide))
			print "hi"
			#self.peptides = SubGroups.reduce(self.peptides)
		self.pair_dict = self.make_pairs()
	
	'''Takes a set of peptides and breaks it down into a set of sequential pairs'''
	def make_pairs(self):
		sequential_pairs = defaultdict(int) #so that all values initialize to 0
		
		# Scan through all the peptides in the docuemtns
		for peptide in self.peptides:
			for i in range(0, len(peptide)-1):
				#count the number of times each pair appears in peptide set
				pair = (peptide[i], peptide[i+1])
				if pair in sequential_pairs.keys():
					sequential_pairs[pair] += 1
				else:
					sequential_pairs[pair] = 1
					
				#now add in the first pair
				if peptide[i] in self.valued_spheres.keys():
					self.valued_spheres[peptide[i]] += 1
				else:
					self.valued_spheres[peptide[i]] = 1
				
				#then add the second pair
				if peptide[i+1] in self.valued_spheres.keys():
					self.valued_spheres[peptide[i+1]] += 1
				else:
					self.valued_spheres[peptide[i+1]] = 1
		return sequential_pairs
	
	'''Takes a sequence pair object and maps the surface residues to the number of
	times they appear in the sequential peptide pairs'''
	def value_pairs(self):
		pair_values = defaultdict(int)
		for pair in self.pair_dict.keys():
			if pair in self.surface_residues:
				#for now we are just adding the pairs if multiple contacts come in
				pair_values[pair] += pair_dict[pair]
		return pair_values
		