#!/usr/bin/python
from single_subject_sl import do_searchlight
from mvpa2.cmdline.helpers import arg2ds
from mvpa2.suite import *

filename = sys.argv[1]
radius   = int(sys.argv[2])

# set up classifiers to try out
clfs = {
        'RidgeRegression': RidgeReg(),
        'LinearSVM': LinearNuSVMC(probability=1,
                      enable_ca=['probabilities']),
        'RBFSVM': RbfNuSVMC(probability=1,
                      enable_ca=['probabilities']),
        'SMLR': SMLR(lm=0.01),
        'Logistic Regression': PLR(criterion=0.00001),
        '3-Nearest-Neighbour': kNN(k=3),
        '10-Nearest-Neighbour': kNN(k=10),
        'GNB': GNB(common_variance=True),
        'GNB(common_variance=False)': GNB(common_variance=False),
#        'LDA': LDA(),
#        'QDA': QDA(),
        }
for name,clf in clfs.iteritems():
	output_basename = os.path.join('{}_r{}_c-{}'.format(filename, radius,name))
        print name
	clf.space= 'condition'
        ds = arg2ds([filename])
        do_searchlight(ds, radius, output_basename,False,clf)

