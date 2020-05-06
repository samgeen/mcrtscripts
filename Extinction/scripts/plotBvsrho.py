'''
Do scatter plots of B versus rho
Sam Geen, July 2019
'''

from startup import *

from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import matplotlib.patheffects as pe

def makeBvsrhohistogram(snap):
    amr = snap.amr_source(["rho","B-left","B-right"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    b = 0.5*(cells["B-left"]+cells["B-right"])  
    ud = snap.info["unit_density"].express(C.g_cc)
    ul = snap.info["unit_length"].express(C.cm)
    ut = snap.info["unit_time"].express(C.s)
    unit = np.sqrt(4.0*np.pi*ud*(ul/ut)**2)*1e6 # microGauss
    Bmag = np.sqrt((b**2).sum(axis=1))*unit
    rhos = cells["rho"]
    histogram = np.histogram2d(rho,Bmag,bins=100,weights=rhos)
    return histogram

if __name__=="__main__":
    sim = imfsims[0]
    snap = sim.Snapshots()[10]
    hist = makeBvsrhohistogram(snap)