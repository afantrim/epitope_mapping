#!/Users/student/anaconda/bin/python
'''
.. module:: PDBUtil
   :platform: Unix, Windows
   :synopsis: Contains miscellaneous utilities for dealing with parsed PDB files. *Consider putting PDBParse in here.*

.. moduleauthor:: Amelia F. Antrim <amelia.f.antrim@gmail.com>

'''

import numpy
from Bio.PDB import *
from collections import defaultdict

class PDBUtil:

    @classmethod
    def calc_residue_dist(cls, residue1, residue2):
    	'''
    	Calculates distance between two residues. Author: Peter Cock, University of Warwick.

    	:param residue1: First residue to calculate distance between.
    	:type residue1: `Bio.PDB.Residue <http://biopython.org/DIST/docs/api/Bio.PDB.Residue-pysrc.html>`_ object.
    	:param residue2: Second residue to calculate distance between.
    	:type residue2: `Bio.PDB.Residue <http://biopython.org/DIST/docs/api/Bio.PDB.Residue-pysrc.html>`_ object.
   
    	'''
    	if "CB" in residue1:
    		atom1 = "CB"
    	else: 
    		atom1 = "CA"
    	if "CB" in residue2:
    		atom2 = "CB"
    	else:
    		atom2 = "CA"
        diff_vector = residue1[atom1].coord - residue2[atom2].coord
        return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

    @classmethod
    def calc_dist_matrix(cls, chain_one, chain_two):
    	'''
    	Creates a nested dictionary/tree containing the distances between all residues in 2 separate
    	chains.

    	:param chain_one: First residue to calculate distance between.
    	:type chain_one: `Bio.PDB.Chain <http://biopython.org/DIST/docs/api/Bio.PDB.Chain.Chain-class.html>`_ object.
    	:param chain_two: Second residue to calculate distance between.
    	:type chain_two: `Bio.PDB.Chain <http://biopython.org/DIST/docs/api/Bio.PDB.Chain.Chain-class.html>`_ object.
   
    	'''
    	print chain_one
    	print chain_two
        ''' Returns a matrix of C-alpha distances between two chains'''
        #answer = numpy.zeros((len(chain_one),len(chain_two)), numpy.float)
        tree = lambda: defaultdict(tree)
        answer = tree()
        
        # Then use nested dictionaries
        for row, residue_one in enumerate(chain_one):
            for col, residue_two in enumerate(chain_two):
                answer[residue_one.id[1]][residue_two.id[1]] = PDBUtil.calc_residue_dist(residue_one, residue_two)
        return answer