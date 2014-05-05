#!/Users/student/anaconda/bin/python

'''
.. module:: PrecisionRecall
   :platform: Unix, Windows
   :synopsis: Makes a precision recall plot from the true and false positives.

.. moduleauthor:: Amelia F. Antrim <amelia.f.antrim@gmail.com>

'''

import numpy as np
import pylab as pl
from sklearn import metrics

class PrecisionRecall:
    
    '''
    .. function: __init__()
    
    Initialize the dictionaries of true and predicted positives
    '''
    def __init__(self):
        self.y_true = {}
        self.y_pred = {}

    def add_threshhold(self, threshold, estimate, true_pos):
    	'''
    	Adds a new set of values for the estimate and the true positives
    	
    	:param threshold: Value of parameter being examined
    	:param estimate: Estimated true positives for a given threshold
    	:param true_pos: True positives for a given threshold
    	'''
        self.y_true[threshold] = true_pos
        self.y_pred[threshold] = estimate
    
    # Draws the precision-recall plot
    def make_plot(self, precision, recall):
    	'''
    	Adds a new set of values for the estimate and the true positives
    	
    	:param precision: Precision values.
    	:type precision: float[]
    	:param recall: Recall values.
    	:type recall: float[]
    	'''
        # Plot the arrays of true and predicted positives on the precision-recall plots
        # Calculate the area under the curve
        pl.clf()
        pl.plot(recall, precision, label='Precision-Recall curve')
        #area = metrics.auc(precision, recall)
        #print("Area Under Curve: %0.2f" % area)
        pl.xlabel('Recall')
        pl.ylabel('Precision')
        pl.ylim([0.0, np.max(precision)+0.05])
        pl.xlim([0.0, np.max(recall)+0.05])
    	#pl.title('Precision-Recall example: AUC=%0.2f' % area)
        pl.legend(loc="lower left")
        pl.show()
        
    def make_roc(self, precision, recall):
    	pass