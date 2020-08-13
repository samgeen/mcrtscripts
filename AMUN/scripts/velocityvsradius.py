'''
Plot velocity vs radius for different wind bubbles
Sam Geen, June 2020
'''

from startup import *

import findproperties, starrelations

from scipy.interpolate import interp1d

maxwindradiusatstarpos = Hamu.Algorithm(findproperties.maxwindradiusatstarpos)

def runforsim(simname):
    print "Running for simulation", simname
    nprocs = 1
    sim = hamusims[simname]
    snap = sim.Snapshots()[0]
    boxlen = snap.RawData().info["boxlen"]
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
    medt = 0.5*(newt[1:nsample] + newt[0:nsample-1])
    # Subtract creation time of star
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    medt -= tcreated
    return medr, vs, medt

def plot(simnames):
    # Get data
    results = {}
    for simname in simnames:
        results[simname] = runforsim(simname)
    # Observed data
    pabstr = 2.0 # pc
    pabstv = 13.0 # km/s
    pabsttime = 0.2 # Myr
    # Plot velocity vs time
    plt.clf()
    for simname in simnames:
        r, v, t = results[simname]
        line = "-"
        if "DENSE" in simname:
            line = "--"
        plt.plot(r,v,color=linestyles.Colour(simname),linestyle=line,label=linestyles.Label(simname))
    plt.plot(pabstr,pabstv,"ok",label="Pabst+ (2019)")
    plt.legend(frameon=False,fontsize="x-small")
    plt.xlabel("Maximum Wind Radius / pc")
    plt.ylabel("Velocity / km/s")
    plt.savefig("../plots/velocityvsmaxradius.pdf")
    # Make radius vs time plot
    plt.clf()
    for simname in simnames:
        r, v, t = results[simname]
        line = "-"
        if "DENSE" in simname:
            line = "--"
        plt.plot(t[t >= 0],r[t >= 0],color=linestyles.Colour(simname),linestyle=line,label=linestyles.Label(simname))
    plt.plot(pabsttime,pabstr,"ok",label="Pabst+ (2019)")
    plt.legend(frameon=False,fontsize="x-small")
    plt.ylabel("Maximum Wind Radius / pc")
    plt.xlabel("Age of Star / Myr")
    plt.savefig("../plots/maxradiusvstime.pdf")

if __name__=="__main__":
    simnames = ["UVWINDPRESS_30"]#,"UVWINDPRESS_60","UVWINDPRESS_120","UVWIND_120_DENSE"]
    plot(simnames)
