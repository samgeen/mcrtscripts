'''
Pymses 4.0 is very inefficient at giving the time, so put it in Hamu
Sam Geen, October 2016
'''

import Hamu

from pymses.utils import constants as C

def _Myr(snap):
    '''
    Returns snapshot time (in Myr)
    '''
    t = snap.info["time"] * snap.info["unit_time"].express(C.Myr)
    return t

# TODO: Include other units as a wrapper around this function
Myr = Hamu.Algorithm(_Myr)
