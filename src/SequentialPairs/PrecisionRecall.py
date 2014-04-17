#!/Users/ameliaantrim/anaconda/bin/python

'''
Makes a precision recall plot from the true and false positives.

NOTE: So far in this context, I use this module to assess the threshhold of what
number of contact pairs are necessary to consider a certain residue a "positive".
'''

import numpy as np
from sklearn.metrics import precision_recall_curve

class PrecisionRecall:
    
    # Initialize the arrays of true and predicted positives
    def __init__(self):
        self.y_true = []
        self.y_pred = []

    def add_threshhold(self, threshhold, true_pos):
        # y_true will hold the numbers of true positives predicted for each
        # y_true will hold the numbers of 
        self.y_true.append(len(true_pos))
        self.y_true.append(len(estimate))
        
    def make_plot(self):
        # Plot the arrays of true and predicted positives on the precision-recall plots
        precision, recall, threshholds = precision_recall_curve(y_true, y_scores)