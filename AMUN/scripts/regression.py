'''
Write the correlation coefficient to file for the scatter plots
Sam Geen, February 2018
'''

from startup import *

def writecoeff(x,y,filename):
    r = np.corrcoef(x,y)
    rxy = r[0,1] # Symmetric matrix where on-axis values are 1
    print("Writing rxy=",str(rxy),"to file",filename)
    f = open(filename,"w")
    f.write(str(rxy))
    f.close()
