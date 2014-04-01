#!/Users/ameliaantrim/anaconda/bin/python

import sys
import math

'''Parses a PDB file. 
Contents:
	self.file = the PDB file containing the protein information
	self.Length = the number of amino acids in the chain
	self.seq = amino acid sequence of chain
	self.StartingResidue = first residue in chain
	self.ContactMap = contact residues, defined as the distance between
		them being no more than "cuttoff" angstroms (here, 8)
	self.DistanceMap = the distances between contact pair residues
	self.Pairs = contact pairs from the contact map and structure info
	self.Distance = distances between the above contact pairs

Warning: Only works on a single chain of the protein. A discontinuous protein warning
often suggests that the protein in question has multiple chains. 
'''
class PDBParse:
    def __init__(self, File):
    	self.file = File
        self.Coord = []
        ResNumList =[]
        AAConverter = {"ALA": "A","GLU": "E","GLN": "Q","ASP": "D","ASN": "N","LEU": "L",\
        "GLY": "G","LYS": "K","SER": "S","VAL": "V","ARG": "R","THR": "T","PRO": "P",\
        "ILE": "I","MET": "M","PHE": "F","TYR": "Y","CYS": "C","CYX": "C","CYD": "C",\
        "TRP": "W","HID": "H","HIE": "H","HIS": "H","HIP":"H"}
        self.seq = ""
        for i in open(File).readlines():
            if i[:4] == "ATOM":
                l = LineParse(i)
                ResCheck = False
                if l.ResName == "GLY":
                    if l.AtomName == "CA":
                        self.Coord.append(l)
                        self.seq += AAConverter[l.ResName]
                        ResNumList.append(l.ResNum)
                        ResCheck = True
                else:
                    if l.AtomName == "CB":
                        self.Coord.append(l)
                        self.seq += AAConverter[l.ResName]
                        ResNumList.append(l.ResNum)
                        ResCheck = True
                if len(ResNumList) > 1 and ResCheck:
                    if ResNumList[-1]-ResNumList[-2] != 1:
                        sys.stderr.write("This protein is not continuous. Check at \
                        %d\n"%ResNumList[-1])
                        sys.exit()
        self.Length = len(self.seq)
        self.StartingResidue = ResNumList[0]

        cutoff = 8
        self.ContactMap = []
        self.DistanceMap = []
        for i in xrange(self.Length):
            if self.seq[i] == "P":
                CM_line = [0]*self.Length
                DM_line = [0]*self.Length
            else:
                CM_line = []
                DM_line = []
                for j in xrange(self.Length):
                    if self.seq[j] == "P":
                        CM_line.append(0)
                        DM_line.append(0)
                    else:
                        a = self.Coord[i]; b= self.Coord[j]
                        distance = math.sqrt((a.X-b.X)**2+(a.Y-b.Y)**2+(a.Z-b.Z)**2)
                        if distance < cutoff: CM_line.append(1)
                        else: CM_line.append(0)
                        DM_line.append(distance)
            self.ContactMap.append(CM_line)
            self.DistanceMap.append(DM_line)

		#get contact pairs within molecule 
        self.Pairs = []
        for i in xrange(self.Length-1):
            for j in xrange(i,self.Length):
                if i==j: pass
                else:
                    if self.ContactMap[i][j]: self.Pairs.append((i,j))
		
		#get distances between contact pairs within molecule
        self.Distances = {}
        for i in xrange(self.Length-1):
            for j in xrange(i,self.Length):
                if i==j: pass
                else:
                    self.Distances[(i,j)] = self.DistanceMap[i][j]

'''Parse and represent a single line in a PDB file'''
class LineParse:
    def __init__(self, line):
        self.AtomName = line[12:16].strip()
        self.ResName = line[17:20].strip()
        self.ResNum = int(line[22:26])
        self.X = float(line[30:38])
        self.Y = float(line[38:46])
        self.Z = float(line[46:54])


def main():
	pdb_file = sys.argv[1]
	pdb = PDBParse(pdb_file)
	print pdb.Pairs
	#print pdb.Distances
	#print pdb.ContactMap
	print pdb.seq
	print pdb.Coord[0]

if __name__ == "__main__":
	main()