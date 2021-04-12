'''
Pymses 4.0 is very inefficient at giving the time, so put it in Hamu
Sam Geen, October 2016
'''

import HamuLite as Hamu

nopymses = False
try:
    from pymses.utils import constants as C
except:
    nopymses = True
    print("No pymses available, skipping Myr conversion")

def _Myr(snap):
    '''
    Returns snapshot time (in Myr)
    '''
    if nopymses:
        t = snap.Time()
    else:
        t = snap.info["time"] * snap.info["unit_time"].express(C.Myr)
    return t

# TODO: Include other units as a wrapper around this function
Myr = Hamu.Algorithm(_Myr)
