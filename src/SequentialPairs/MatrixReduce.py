#!/Users/student/anaconda/bin/python

'''
.. module:: MatrixReduce
   :platform: Unix, Windows
   :synopsis: Look up 2 pairs of amino acids according to a dictionary that can map them to a certain score based on how close the pairs are.

.. moduleauthor:: Amelia F. Antrim <amelia.f.antrim@gmail.com>

'''

import Blosum

class MatrixReduce:
	
	@classmethod
	def score_pair(cls, true_pair, pair, scoring_matrix):
		'''
		Given two pairs of residues, calculates a substitutability
		index for these pairs given a scoring matrix.

    	:param true_pair: The pair of residues in the true protein.
    	:type true_pair: Tuple of `Bio.PDB.Residue <http://biopython.org/DIST/docs/api/Bio.PDB.Residue-pysrc.html>`_ objects.
    	:param pair: Pair of residues for which we assess substitutability.
    	:type pair: Tuple of `Bio.PDB.Residue <http://biopython.org/DIST/docs/api/Bio.PDB.Residue-pysrc.html>`_ objects.
   		:param scoring_matrix: Matrix containing substitutability scores.
   		:type scoring_matrix: Nested dictionary mapping 2 residues to their substitutability score.
		:returns: float-- the substitutability score of the pairs 
		
		'''
	
		if scoring_matrix == None:
			scoring_matrix = Blosum.CSERZO
			
		# Check both amino acids in the pair against the true pair
		if pair[0].get_resname() in scoring_matrix[true_pair[0]].keys():
			score = scoring_matrix[true_pair[0]][pair[0].get_resname()]
		else:
			score = scoring_matrix[pair[0].get_resname()][true_pair[0]]
		if pair[1].get_resname() in scoring_matrix[true_pair[1]].keys():
			score += scoring_matrix[true_pair[1]][pair[1].get_resname()]
		else:
			score += scoring_matrix[pair[1].get_resname()][true_pair[1]]
		
		# Count the score if it is high enough to be a possible true positive
		return score