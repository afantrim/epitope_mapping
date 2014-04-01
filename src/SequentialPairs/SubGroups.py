#!/Users/ameliaantrim/anaconda/bin/python

'''Dictionary that converts amino acids into the function groups defined by Mapitope'''


'''Represents a string of amino acids (a peptide or protein) in both original form and
reduced-alphabet form'''
class SubGroups:
	'''def __init__(self, aa_string):
		self.aa_string = aa_string
		self.reduced_string = self.reduce()'''
	
	'''Takes a string of amino acids and maps it to the reduced amino acid alphabet
	published in Mapitope'''
	@staticmethod
	def reduce(aa_string):
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
		for i in range(0, len(aa_string)):
			aa_string[i] = aaDict[aa_string[i]]
		return aa_string
	
	@staticmethod
	def residue_reduce(generator):
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
			print aa_string
			print new_str
		return new_str