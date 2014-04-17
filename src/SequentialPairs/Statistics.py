#!/Users/ameliaantrim/anaconda/bin/python

import sys
import os

class Statistics:
    def __init__(self, predicted_contacts, true_contacts):
        self.positives(predicted_contacts, true_contacts)

    def positives(self, predicted_contacts, true_contacts):
        # Make a hash table to check the validity of positives
        positives = {} # Doesn't need to be a self b/c only used internally for validation
        self.true_positives = []
        self.false_positives = []
        self.true_negatives = []

        # Add all the true contact pairs to the hash table
        # Using a 0/1 flag to indicate whether the residue was found or not
        for residue in true_contacts:
            positives[residue] = 0
        
        # Iterate through the predicted positives and check them against the true values
        for residue in predicted_contacts:
            # If the predicted positive is in the true positive hash
            if residue in positives.keys():
                # Mark it as estimated 
                self.true_positives.append(residue)
                # Change the flag to "discovered
                positives[residue] += 1
            # Else it is a false positive
            else:
                self.false_positives.append(residue)

        # Now check whether the true positives were predicted or not
        for residue in positives.iterkeys():
            if positives[residue] == 0:
                self.false_negatives.append(residue)

        # If we want to find true negatives, we would iterate through the entire structure
        # and record everything that wasn't already in any of these categories
        # Check in about whether this is necessary11