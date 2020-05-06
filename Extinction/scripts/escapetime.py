'''
Time taken in each simulation from first star to photons escaping cloud
Sam Geen, March 2018
'''

from startup import *

import scipy.interpolate

import plotproperties, starrelations

def tescape(simname,rescape=10.0):
    '''
    Find escape time for a simulation
    rescape = radius to use for escape / pc
    '''
    t,r = plotproperties.radius(simname)
    fr = scipy.interpolate.interp1d(r,t,kind="linear")
    tesc = fr(rescape)
    return tesc

def runforsim(simname):
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    tesc = tescape(simname,rescape=5.0)
    return tesc - tcreated

def run(simnames):
    tescs = []
    for simname in simnames:
        tescs.append(runforsim(simname))
    print tescs

if __name__=="__main__":
    run(imfsims)
