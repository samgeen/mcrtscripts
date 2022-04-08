'''
Plot the solid angle subtended by the wind bubble at each radius
Sam Geen, April 2022
'''

from startup import *

import rayprof

from collections import OrderedDict

nprofs = rayprof.nprofs
nray = rayprof.nray

def findstarpos(snap,useMyr=False):
    stellar = stellars.FindStellar(snap)
    if len(stellar.mass) == 0:
        print("No stellar objects in snapshot")
        return None
    imax = np.where(stellar.mass == np.max(stellar.mass))[0][0]
    sinkid = stellar.sinkid[imax]-1
    sink = sinks.FindSinks(snap)
    boxlen = snap.RawData().info["boxlen"]
    pos = np.array([sink.x[sinkid],sink.y[sinkid],sink.z[sinkid]])/boxlen
    return pos


def plotsolidangles(simname,times):
    '''
    Like prof but with a higher gradient
    '''
    global nray, nprofs,rays,rads
    powfile = None
    xlims = None
    #sim = Hamu.Simulation(simname)
    sim = hamusims[simname]
    # HACK
    print("Running for sim", simname)
    snaps = []
    for timetuple in times:
        snaps.append(timefuncs.findsnapshot(sim,timetuple))
    plt.clf()
    lines = ["-","--",":"]
    iline = -1
    # Find ray where shock is most extended
    for timetuple, snap in zip(times,snaps):
        iline += 1
        line = lines[iline]
        time, tunit = timetuple
        tstr = str(time)
        if "Myr" in tunit:
            tstr = str(time)+" Myr"
        starpos = findstarpos(snap)
        if starpos is None:
            continue # skip this iteration
        centre = starpos
        #r,nHs   = rayprof.makeprofsHamu(snap,"nH",centre)
        r,Ts    = rayprof.makeprofsHamu(snap,"T",centre)
        r,xHIIs = rayprof.makeprofsHamu(snap,"xHII",centre)
        radii = r[0:nray]
        #nHs = np.reshape(nHs,(nprofs, nray)).T
        Ts = np.reshape(Ts,(nprofs, nray)).T
        xHIIs = np.reshape(xHIIs,(nprofs, nray)).T
        iraymax = None
        maxr = 0.0
        Tthresh = 1e6
        xHIIthresh = 0.1
        windinray = radii*0.0
        colours = {True: "k", False:"r"}
        typelines = []
        for doHII in [True, False]:
            colour = colours[doHII]
            for iray in range(0,nray):
                raymask = np.where(Ts[:,iray] > Tthresh)
                if doHII:
                    raymask = np.where((Ts[:,iray] > Tthresh) | (xHIIs[:,iray] > xHIIthresh))
                windinray[raymask] += 1
            windinray *= 4.0 * np.pi / nray
            mask = np.where(windinray > 0.0)
            radii = radii[mask]
            windinray = windinray[mask]
            l1, = plt.plot(radii,windinray,linestyle=line,label=tstr,color=colour)
            tstr = None
            typelines.append(l1)
    plt.xlabel("Radius / pc")
    sub = "wind"
    if doHII:
        sub = "HII"
    plt.ylabel("$\Omega_{\mathrm{wind,HII}}(R)$")
    leg1 = plt.legend(frameon=False)
    leg2 = plt.legend(typelines,labels=["HII","wind"],frameon=False,loc="center right")
    plt.gca().add_artist(leg1)
    suffix = ""
    plt.xscale("linear")
    plt.yscale("linear")
    yticks = [i*np.pi for i in [0,1,2,3,4]]
    yticklabels = ["0","$\pi$","$2 \pi$","$3 \pi$","$4 \pi$"]
    plt.gca().set_yticks(yticks)
    plt.gca().set_yticklabels(yticklabels)
    #plt.ylim([1e-14,1e10])
    filename = "../plots/windandHIIsolidangle_"+simname+".pdf"
    print("Writing "+filename)
    plt.savefig(filename)
    

if __name__=="__main__":
    timetuples = [(t,"MyrFirstStar") for t in (0.1,0.2,0.3)]
    simname = "SEED1_35MSUN_CDMASK_WINDUV"
    plotsolidangles(simname,timetuples)
