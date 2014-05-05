#!/Users/student/anaconda/bin/python
'''
.. module:: Statistics
   :platform: Unix, Windows
   :synopsis: Generates statistics about results. 

.. moduleauthor:: Amelia F. Antrim <amelia.f.antrim@gmail.com>

'''


import sys
import os

class Statistics:
    '''
    Generates lists of true and false positives/negatives given predictions and true answers.
    
    :param predicted_contacts: estimated positives
    :type predicted_contacts: array
    :param true_contacts: list of true positives (whether predicted or not)
    :type true_contacts: array
    :param antigen_chain: antigen structure
    :type antigen_chain: Bio.PDB.Chain object. 
    
    '''
    def __init__(self, predicted_contacts, true_contacts, antigen_chain):
        self.positives(predicted_contacts, true_contacts, antigen_chain)

    def positives(self, predicted_contacts, true_contacts, antigen_chain):
        # Make a hash table to check the validity of positives
        positives = {} # Doesn't need to be a self b/c only used internally for validation
        self.true_positives = []
        self.false_positives = []
        self.false_negatives = []
        
        print "true contacts"
        print true_contacts
        # Add all the true contact pairs to the hash table
        # Using a 0/1 flag to indicate whether the residue was found or not
        for residue in true_contacts:
            positives[residue] = 0
        
        # Iterate through the predicted positives and check them against the true values
        for residue in predicted_contacts:
            # If the predicted positive is in the true positive hash
            if residue in positives.keys():
                # Mark it as estimated 
                #print "\nTRUE: Comparing " + str(residue) + " and " + str(positives[residue])
                self.true_positives.append(residue)
                # Change the flag to "discovered
                positives[residue] += 1
            # Else it is a false positive
            else:
                self.false_positives.append(residue)
                #print "\nFALSE: Comparing " + str(residue)

        # Now check whether the true positives were predicted or not
        for residue in positives.iterkeys():
            if positives[residue] == 0:
                self.false_negatives.append(residue)
        
        return self.true_positives, self.false_positives

        # If we want to find true negatives, we would iterate through the entire structure
        # and record everything that wasn't already in any of these categories
        # Check in about whether this is necessary11