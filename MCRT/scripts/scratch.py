'''
TESTBED FOR STUFF
'''

import customplot
import matplotlib.pyplot as plt

bee = "bee"
plt.plot([1,2],[3,4],label="bah"+bee)
leg = plt.legend(ncol=2)
for text in leg.get_texts():
    print text
    text = text+text
plt.savefig("scratch.pdf")
