#!/Users/ameliaantrim/anaconda/bin/python

'''Dictionary that converts amino acids into the function groups defined by Mapitope'''
import Bio.PDB

'''Represents a string of amino acids (a peptide or protein) in both original form and
reduced-alphabet form'''
class SubGroups:
	
	'''
	Takes  a residue object and changes its resname to the mapping of that
	residue onto the reduced alphabet.
	'''
	@staticmethod
	def residue_reduce(residue):
		aaDict = {
			'R': 'B',
			'K': 'B',
			'E': 'J',
			'D': 'J',
			'S': 'O',
			'T': 'O',
			'L': 'U',
			'V': 'U',
			'I': 'U',
			'Q': 'X',
			'N': 'X',
			'W': 'Z',
			'F': 'Z',
			'A': 'A',
			'C': 'C',
			'G': 'G',
			'H': 'H',
			'M': 'M',
			'P': 'P',
			'Y': 'Y'
		}
		#print type(residue)
		# Make sure we are working with a residue
		if type(residue) is not Bio.PDB.Residue.Residue:
		    print "Incorrect object passed in ot residue_reduce."
		    exit(1)
		
		# Convert the resname in the residue object to its reduced name
		residue.resname = aaDict[residue.get_resname()]
		return residue.resname
	
	'''
    Takes in an amino acid sequence as a string and returns the corresponding string
    in the reduced alphabet. 
	'''
	@staticmethod
	def str_reduce(aa_string):
		aaDict = {
			'R': 'B',
			'K': 'B',
			'E': 'J',
			'D': 'J',
			'S': 'O',
			'T': 'O',
			'L': 'U',
			'V': 'U',
			'I': 'U',
			'Q': 'X',
			'N': 'X',
			'W': 'Z',
			'F': 'Z',
			'A': 'A',
			'C': 'C',
			'G': 'G',
			'H': 'H',
			'M': 'M',
			'P': 'P',
			'Y': 'Y'
		}
		new_str = ""
		for i in range(0, len(aa_string)):
			if aa_string[i] in aaDict.keys():
				new_str += aaDict[aa_string[i]]
			else:
				new_str += aa_string[i]
		return new_str