'''
Generic function that returns time series for an input function and arguments
Sam Geen, December 2017
'''

import numpy as np
import snaptime

import multiprocessing as mp

ERRCODE = None
IGNOREERRORS = False

class ErrFunc(object):
    # Wrap a function around a possible error
    def __init__(self,func):
        self._func = func

    def __call__(self,snap,*args,**kwargs):
        try:
            ret = self._func(snap,*args,**kwargs)
        except ValueError:
            print("ERROR FOUND, IGNORING OUTPUT...")
            return ERRCODE

def timefunc(sim,func,noarray=False,verbose=False,processes=1,*args,**kwargs):
    times = []
    vals = []
    errfunc = ErrFunc(func)
    torun = func
    if IGNOREERRORS:
        torun = errfunc
    if processes == 1:
        if verbose:
            print(sim.Name())
        for snap in sim.Snapshots():
            t = snaptime.Myr(snap)
            v = torun(snap,*args,**kwargs)
            if verbose:
                print("t,v:", t, v)
            times.append(t)
            vals.append(v)
    else:
        pool = mp.Pool(processes=processes)
        times = pool.map(snaptime.Myr, sim.Snapshots())
        vals  = pool.map(torun,        sim.Snapshots())

    times = np.array(times)
    if not noarray:
        vals = np.array(vals)
    if IGNOREERRORS:
        mask = np.where(vals is not ERRCODE)
        times = times[mask]
        vals = vals[mask]
    if verbose:
        print(times,vals)
    return times, vals

def findsnapshot(sim,timetuple):
    # Function for finding at a time with different criteria (e.g. Myr, since first star)
    simname = sim.Name()
    # Time stuff
    try:
        snap = sim.Snapshots()[0]
    except:
        print("ERROR FINDING ANY SNAPSHOTS IN", simname)
        print(sim.Folder())
        raise ValueError
    myr   = snap.RawData().info["unit_time"].express(C.Myr)
    time, timeunits = timetuple
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    tcreatedcode = tcreated/myr
    if timeunits == "Myr":
        time /= myr
    if timeunits == "MyrFirstStar":
        time += tcreated # start from time first star created
        time /= myr
    if timeunits == "code":
        pass # already ok
    if timeunits == "codeFirstStar":
        time += tcreated
    if timeunits == "outputNumber":
        outsnaps = {snap.OutputNumber():snap for snap in sim.Snapshots()}
        time = outsnaps[int(time)].Time()
    return sim.FindAtTime(time)

