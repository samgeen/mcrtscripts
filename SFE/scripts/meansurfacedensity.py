'''
Find the mean surface density in a cloud
Sam Geen, November 2016
'''

from startup import *
import columndensity, freefalltime, linestyles
import Hamu

def run(sim,los):
    # Find snapshot at tff (in code units)
    tff = freefalltime.Tff(sim)
    snap = sim.FindAtTime(tff)
    surf = runinsnap(snap,los)
    print surf,"/ Msun/pc^2 in "+sim.Name() 
    return surf

def _meansurfacedensityinsnap(snap,los):
    # Get density map
    nhlow = AktoNH(0.1)
    dmap = columndensity.DensityMap(snap,los,NHlow=nhlow)
    # Find properties of map
    mass = np.sum(dmap.Mass())
    area = dmap.Area()
    return mass / area
runinsnap = Hamu.Algorithm(_meansurfacedensityinsnap)

def runall():
    sims = [Hamu.Simulation(x+"-NRT") for x in ["L","M","S","XS"]]
    results = {}
    for sim in sims:
        vals = []
        for los in ["x","y","z"]:
            vals.append("%.1f" % run(sim,los))
        vals.sort()
        results[sim.Name()] = vals[1]+"^{"+vals[2]+"}_{"+vals[0]+"}"
    keys = results.keys()
    keys.sort()
    for key in keys:
        print key, results[key]

if __name__=="__main__":
    runall()

