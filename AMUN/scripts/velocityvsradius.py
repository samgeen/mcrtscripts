'''
Plot velocity vs radius for different wind bubbles
Sam Geen, June 2020
'''

from startup import *

from scipy.interpolate import interp1d

maxwindradiusatstarpos = Hamu.Algorithm(findproperties.maxwindradiusatstarpos)

def runforsim(simname):
    print "Running for simulation", simname
    nprocs = 10
    sim = hamusims[simname]
    # Get max wind radii in each snapshot
    t,r = timefuncs.timefunc(sim,maxwindradiusatstarpos,verbose=True,processes=nprocs)
    # Resample at uniform times to get velocity estimates
    rfunc = interp1d(t,r)
    nsample = 1000
    newt = np.linspace(t[0],t[-1],nsample)
    newr = rfunc(newt)
    dt = newt[1] - newt[0]
    vs = (newr[1:nsample] - newr[0:nsample-1]) / dt
    # Get median radius in the sample
    medr = 0.5*(newr[1:nsample] + newr[0:nsample-1])
    return medr, vs

if __name__=="__main__":
    runforsim("UVWINDPRESS_30")