'''
Calculate the exinction of sinks in the cloud
Sam Geen, June 2016
'''

from startup import *

import hydrofuncs
import columndensity, sinks

lostonum = {'x':0,'y':1,'z':2}

def SinkExtinction(snap,los):
    '''
    For each sink, returns two arrays, showing extinction viewed 
      from either face
    '''
    up = columndensity.ups[los]
    across = columndensity.acrosses[los]
    sinklist = sinks.FindSinks(snap.hamusnap)
    if sinklist.IsEmpty():
        return np.array([]), np.array([])
    nray = 512
    amr = hydrofuncs.amr_source(snap,"rho")
    boxlen = snap.info["boxlen"]
    ray = np.zeros((nray,3)) # We overwrite this each time
    def _FindExtinction(snap,x,y,z,los,frontToBack=True):
        x /= boxlen
        y /= boxlen
        z /= boxlen
        ray[:,0] = x
        ray[:,1] = y
        ray[:,2] = z
        pos = [x,y,z]
        losstart = pos[lostonum[los]]
        # What is the end point?
        losend = 0.0
        if frontToBack:
            losend = 1.0
        dlos = np.abs(losstart-losend)/float(nray)
        # Make sure the start is smaller than the end
        if losstart > losend:
            losstart, losend = losend, losstart
        ray[:,lostonum[los]] = np.arange(losstart,losend,dlos)
        # Now sample along ray
        dset = pymses.analysis.sample_points(amr,ray)
        rho = hydrofuncs.scale_by_units(snap,"rho")(dset)
        lrayinpc = np.abs(losstart-losend)*boxlen*pcincm/float(nray)
        NHout = np.sum(rho) * lrayinpc
        Akout = NHtoAk(NHout)
        return Akout
    FindExtinctionHamu = Hamu.Algorithm(_FindExtinction)
    # Run through sinks
    lsinks = sinklist.Length()
    Akftb = np.zeros(lsinks)
    Akbtf = np.zeros(lsinks)
    for i, sx, sy, sz in zip(range(0,lsinks),
                             sinklist.x, sinklist.y, sinklist.z):
        Akftb[i] = FindExtinctionHamu(snap.hamusnap,sx,sy,sz,los)
        Akbtf[i] = FindExtinctionHamu(snap.hamusnap,sx,sy,sz,los,
                                   frontToBack=False)
    return Akftb, Akbtf

SinkExtinctionHamu = Hamu.Algorithm(SinkExtinction)


if __name__=="__main__":
    # Test this
    from pymses.utils import constants as C
    sim = Hamu.Simulation("M-RT")
    for snap in sim.Snapshots():
        utime = snap.RawData().info["unit_time"].express(C.Myr)
        time = snap.Time()*utime
        ftb, btf = SinkExtinctionHamu(snap,'z')
        print ftb, btf, time
