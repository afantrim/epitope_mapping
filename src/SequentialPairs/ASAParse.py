#!/Users/student/anaconda/bin/python

import sys
import math
import os

'''Parse an ASA file to extract contact pairs on the surface of the molecule.'''
class ASAParse:
    def __init__(self, File):
        self.SA = []
        self.SA_dict = {}
        try:
            for line in open(File, 'r'):
                self.SA.append(float(line.split(' ')[1]))
                self.SA_dict[line.split(' ')[0][1::]] = float(line.split(' ')[1])
            print 'Reading solvent accessibility file ... done'
        except:
            print 'Missing solvent accessibility info file.'
            print 'Please provide a two column file that has wildtype ASA info.'
            print 'Use ASAView server http://www.abren.net/asaview/ to get this info \
            from your PDB.'