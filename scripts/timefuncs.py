'''
Generic function that returns time series for an input function and arguments
Sam Geen, December 2017
'''

import numpy as np
import snaptime

def timefunc(sim,func,noarray=False,*args,**kwargs):
    times = []
    vals = []
    for snap in sim.Snapshots():
        times.append(snaptime.Myr(snap))
        vals.append(func(snap,*args,**kwargs))
    times = np.array(times)
    if not noarray:
        vals = np.array(vals)
    return times, vals
