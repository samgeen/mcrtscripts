'''
Find the YSOs in a set of sink particles based on the accretion history
Sam Geen, October 2016
'''

from startup import *

import columndensity, images, sinks

sim = None

def _SinksAboveThresh(snap,los,Aklim=0.1):
    xl = columndensity.IMSIZE
    yl = columndensity.IMSIZE
    dmap = columndensity.DensityMap(snap,los,NHlow=AktoNH(Aklim))
    boxlen = snap.info["boxlen"]
    sinkx, sinky, sinkm = images.ProjectSinks(snap.hamusnap,los)
    def NHAtPos(x,y):
        # Find the column density in im at x,y
        # x,y position is in coordinates (0,boxlen)
        px = (xl*x/boxlen).astype(int)
        py = (yl*y/boxlen).astype(int)
        # Reset minimum density
        dmap.NHlow(AktoNH(0.0))
        return dmap.NH()[px,py]
    # Get list of column densities
    NHs = NHAtPos(sinkx, sinky)
    # Find mass of sinks above NH threshold and below extinction limit
    NHthresh = AktoNH(Aklim)
    return NHs >= NHthresh
SinksAboveThresh = Hamu.Algorithm(_SinksAboveThresh)

def FindYSOMass(times,masses,ysoage,time,mass):
    '''
    Find the mass in YSOs
    times in Myr, stars as cumulative mass in Msun, ysoage in Myr
    '''
    ysomass = mass
    if time > ysoage:
        ysomass = mass - np.interp(t-ysoage,times,stars)
    return ysomass

def _MassInSnap(snap,los,ysoages,Aklim=0.1,persink=False):
    '''
    Find the mass of YSOs in a snapshot
    OBLIGATORY NAMELIST VARIABLE sim: The simulation (to build the sink history)
    ysoages: Age of YSOs in Myr [ array of ages ]
    Aklim: extinction limit to use
    '''
    global sim
    if sim is None:
        print "MUST SET sim VARIABLE ON A NAMELIST LEVEL:", sim
        raise ValueError
    time = snaptime.Myr(snap.hamusnap)
    # Build list of snapshots to compare
    snaps = []
    times = []
    for s in sim.Snapshots():
        stime = snaptime.Myr(s)
        if stime <= time:
            snaps.append(s)
            times.append(stime)
    times = np.array(times)
    # Find all sinks in the snapshot
    currsinks = sinks.FindSinks(snap.hamusnap)
    currids = currsinks.id
    # No sinks, just return
    if len(currids) == 0:
        return ysoages*0.0
    # Find IDs in the Ak limit
    threshids = currids[SinksAboveThresh(snap.hamusnap,los,Aklim)]
    # Make table of sink masses
    masses = np.zeros((currids.max()+1,len(snaps)))
    isnap = -1
    for s in snaps:
        isnap += 1
        sink = sinks.FindSinks(s)
        try:
            a = len(sink.id)
        except:
            import pdb; pdb.set_trace()
        if len(sink.id) > 0:
            masses[sink.id,isnap] += sink.mass
    # Find YSO mass for each sink
    ysomasses = {a: [] for a in ysoages}
    for cid in currids:
        if cid in threshids:
            mass = masses[cid,:]
            for age in ysoages:
                ysomass = mass[-1]
                if time > age:
                    ysomass -= np.interp(time-age,times,mass)
                ysomasses[age].append(ysomass)
    # Done!
    if not persink:
        for age in ysoages:
            ysomasses[age] = np.sum(ysomasses[age])
    return ysomasses
MassInSnap = Hamu.Algorithm(_MassInSnap)
