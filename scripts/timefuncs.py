'''
Generic function that returns time series for an input function and arguments
Sam Geen, December 2017
'''

import numpy as np
import snaptime

import multiprocessing as mp

def timefunc(sim,func,noarray=False,verbose=False,processes=1,*args,**kwargs):
    times = []
    vals = []
    if processes == 1:
        for snap in sim.Snapshots():
            t = snaptime.Myr(snap)
            v = func(snap,*args,**kwargs)
            if verbose:
                print "t,v:", t, v
            times.append(t)
            vals.append(v)
    else:
        pool = mp.Pool(processes=processes)
        ts   = pool.map(snaptime.Myr, sim.Snapshots())
        vals = pool.map(func,         sim.Snapshots())
    times = np.array(times)
    if not noarray:
        vals = np.array(vals)
    return times, vals
