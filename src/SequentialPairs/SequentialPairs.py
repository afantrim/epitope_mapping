#!/Users/student/anaconda/bin/python

'''
.. module:: SequentialPairs
   :platform: Unix, Windows
   :synopsis: Deals with peptide binding datasets and the sequential subsets thereof. 

.. moduleauthor:: Amelia F. Antrim <amelia.f.antrim@gmail.com>

'''

import sys
import math
import os
from collections import defaultdict
from ASAParse import ASAParse
from SubGroups import SubGroups

'''
Represents a set of sequential pairs of amino acids derived from a set of peptide
binding data, maps them to the surface pairs of the antigen
'''
class SequentialPairs:
	
	'''  
    .. function:: __init__(self, peptide_file, reduced)
    
    Represents a set of reduced or non-reduced pairs in the peptide dataset 

    :param peptide_file: File containing list of peptides observed to bind antibody, separated by newline characters.
    :type peptide_file: str (filename)
    :param reduced: Boolean indicating whether or not to reduce the amino acid alphabet.
    :type reduced: boolean or integer (0=False, 1=True)
   
    '''
	def __init__(self, peptide_file, reduced):
		self.peptide_file = peptide_file
		with open(peptide_file, "r") as file:
			self.peptides = file.read().split('\n')
		
		# Convert the amino acid alphabet to the reduced alphabet if desired
		'''if reduced:
			temp = self.peptides
			self.peptides = []
			for peptide in temp:
				self.peptides.append(SubGroups.str_reduce(peptide))'''
		
		# Make pairs
		self.pair_dict = self._make_pairs()
	
	
	def _make_pairs(self):
		'''Takes a set of peptides and breaks it down into a set of sequential pairs.

    	'''
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
		return sequential_pairs
		