'''
Test Bayesian correlation
Sam Geen, March 2018
'''

from startup import *

import syscorr

x = np.random.norm(2,0.1,size=13)
y = np.random.norm(2,0.1,size=13)
chain = np.transpose([x, y])
