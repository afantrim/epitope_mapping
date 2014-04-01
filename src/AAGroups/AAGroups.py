import defaultdict
from __future__ import print_function

#aaDict = defaultdict(char)
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
	#maybe get rid of everything past here and find
	#a way to keep it the same; but this might also
	#be a good check to make sure all AAs are valid
	'A': 'A',
	'C': 'C',
	'G': 'G',
	'H': 'H',
	'M': 'M',
	'P': 'P',
	'Y'; 'Y'
}

#takes newline-delimited file
class AminoAcidSubgroups:
    def __init__(self, file):
    	self.file = file
    	self.reduced_file = make_dict_file
    
    def make_dict_file(self):
        #make new file by appending _grouped to end of PDB name before.txt
        newFile = self.file[0:len(file)-4]+"_grouped.txt"
        try:
            with open(file, 'r') as textFile:
                with open(newFile, 'w') as modFile:
                    for line in textFile: #create simplified peptide for each
                        newString=""
                            for char in line: #translate each character in the peptide
                                newString+=aaDict[char]
                                #write the converted string to the new file
                                modFile.write(newString+'\n')
        except IOError:
            print "The peptide file could not be opened."
        return newFile
            
            
            
            
            
            
            
            
            
            
            
            
            