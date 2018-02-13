'''
Run all plotting tools
Sam Geen, February 2017
'''

from startup import *

import starrelations, correlatestructure

def runall():
    runstarrelations()
    runcorrelatestructure()

def runstarrelations():
    for toplotnow in ["compactness","alltimemax","firstmass","firsttime","nphotons","nphotonstot"]:
        starrelations.run(imfsims,"imf",toplotnow)

def runcorrelatestructure():
    for toplotnow in "F,T,L,M,S,C".split(","):
        for rhofloornow in [10.0,100.0,1000.0]:
            correlatestructure.parameter = toplotnow
            correlatestructure.rhofloor = rhofloornow
            correlatestructure.run(icsims,"ic")

if __name__=="__main__":
    runall()
