import sys, math

# To get ASA data follow these steps:
# 1. Go to http://www.abren.net/asaview/
# 2. Graphical results are displayed with a link to download "relative ASA values (%)" as a text file.
# 3. Parse this text file in ASAParser python code below.


class ASAParse:
    def __init__(self, File):

        self.SA = []
        try:
            for line in open(File, 'r'):
                self.SA.append(float(line.split(' ')[1]))
            #print 'Reading solvent accessibility file ... done'
        except:
            print 'Missing solvent accessibility info file.'
            print 'Please provide a two column file that has wildtype ASA info.'
            print 'Use ASAView server http://www.abren.net/asaview/ to get this info from your PDB.'

        #print self.SA
