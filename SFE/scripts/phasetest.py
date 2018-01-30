'''
Test the gas phase in an output to see how x and T match up
Sam Geen, July 2016
'''

from startup import *

from pymses.filters import CellsToPoints
from pymses.utils import constants as C

def _xvsTinsnap(snap):
    print "Finding x vs T in snap", snap.iout
    amr = snap.amr_source(["rho","P","xHII"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    x = cells["xHII"]
    T = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    return x, T

xvsTinsnap = Hamu.Algorithm(_xvsTinsnap)

def PlotPhase(snap,name):
    print "Plotting x vs T in snap", snap.RawData().iout, "with name", name
    plt.clf()
    x, T = xvsTinsnap(snap)
    # Downsample
    newinds = np.random.choice(len(x),10000)
    x = x[newinds]
    T = T[newinds]
    plt.scatter(x,T)
    plt.xlabel(r"$x_{HII}$")
    plt.ylabel(r"T / K")
    plt.yscale("log")
    plt.savefig("../plots/"+name+"_"+str(snap.RawData().iout).zfill(5)+".pdf")

if __name__=="__main__":
    Hamu.Workspace("HIISFE")
    sim = Hamu.Simulation("L-RT")
    snap = sim.Snapshots()[50]
    PlotPhase(snap,"phase_L-RT")
