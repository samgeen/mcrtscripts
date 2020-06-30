'''
Run all plotting tools
Sam Geen, February 2017
'''

from startup import *

import starrelations, correlatestructure

def runall():
    runcorrelatestructure()
    runstarrelations()

def runstarrelations():
    for toplotnow in ["compactness","alltimemax","firstmass","firsttime","nphotons","nphotonstot","nphotonstff","tracklength"][::-1]:
        starrelations.run(imfsims,"imf",toplotnow)
        starrelations.run(icsims,"ic",toplotnow)

def runcorrelatestructure():
    for toplotnow in "V,F,T,L,M,S,C".split(","):
        for rhofloornow in [10.0,100.0,1000.0]:
            correlatestructure.tff_fact = None
            correlatestructure.parameter = toplotnow
            correlatestructure.rhofloor = rhofloornow
            correlatestructure.run(icsims,"ic")
        for tff_fact in [0.5,1.0,1.5,2.0]:
            correlatestructure.tff_fact = tff_fact
            correlatestructure.parameter = toplotnow
            correlatestructure.rhofloor = 10.0
            correlatestructure.run(icsims,"ic")

if __name__=="__main__":
    runall()
