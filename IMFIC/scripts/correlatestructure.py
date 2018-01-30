'''
Correlate the structure of the cloud to the SFE
Sam Geen, January 2018
'''

from startup import *

import pymses
from pymses.filters import CellsToPoints 
from pymses.utils import constants as C

import triaxfinder, tsfe

#parameter = "F" # Form factor a*b/c2
parameter = "T" # (1-b^2) / (1-g^2)
rhofloor = 100.0

def _findtriax(snap,massweight=True,cutrho=10.0):
    # Find the triaxial structure of the cloud
    # Access cell data
    amr = snap.amr_source(["rho"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    rhos = cells["rho"]*snap.info["unit_density"].express(C.H_cc)
    cellsizes = cells.get_sizes()*snap.info["unit_length"].express(C.pc)
    # Mask out cells below a density threshold (say, 10 atcc)
    mask = rhos > cutrho
    masses = cells["rho"]*cells.get_sizes()**3.0
    masses = masses[mask]
    posns = cells.points[mask,:]
    if not massweight:
        masses *= 0.0
        masses += 1.0
    x = posns[:,0]
    y = posns[:,1]
    z = posns[:,2]
    # Run triaxial finder
    a = np.zeros(1)
    b = np.zeros(1)
    c = np.zeros(1)
    ncells = masses.shape[0]
    a,b,c = triaxfinder.findaxes(masses,x,y,z,ncells)
    return a,b,c
findtriax = Hamu.Algorithm(_findtriax)

def FormParameter(a,b,c):
    '''
    Find a parameter describing the form of the object
    let a,b be most similar and c least similar
    F = a*b/c**2
    F < 1 -> filament
    F > 1 -> disk
    F ~ 1 -> sphere
    '''
    axes = np.sort([a,b,c])
    if parameter == "F":
        # Find value furthest from middle axis
        b = axes[1]
        if np.abs(axes[0]-axes[1]) > np.abs(axes[2]-axes[1]):
            c = axes[0]
            a = axes[2]
        else:
            c = axes[2]
            a = axes[0]
        # Return form parameter
        return a*b/c**2
    if parameter == "T":
        # Kimm & Yi, 2007, eqn 4
        beta = axes[1]/axes[2]
        gamma = axes[0]/axes[2]
        return (1.0-beta**2) / (1.0-gamma**2)

def runforsim(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    # Find form factor at 0.5 t_ff
    tcode = tffcloud_code/2.0
    snap = sim.FindAtTime(tcode)
    a,b,c = findtriax(snap,cutrho=rhofloor)
    F = FormParameter(a,b,c)
    print F
    # Find last TSFE
    snap = sim.Snapshots()[-1]
    sfe = tsfe.tsfeinsnap(snap)
    col = linestyles.Colour(simname)
    return sfe,F,col

def run(simnames,plotname):
    plt.clf()
    tsfes = []
    Fs = []
    cols = []
    for simname in simnames:
        sfe, F, col = runforsim(simname)
        tsfes.append(sfe)
        Fs.append(F)
        cols.append(col)
    tsfes = np.array(tsfes)
    Fs = np.array(Fs)
    plt.scatter(Fs,tsfes*100,s=80,marker="o",c=cols,edgecolors='k')
    if parameter == "F":
        maxF = np.max([np.max(Fs),10.0])
        minF = np.min([np.min(Fs),0.1])
        plt.xscale("log")
        plt.xlabel(r"Filamentary $\leftarrow$ $F$ (\textonehalf $t_{ff}) \rightarrow$ Disky")
    if parameter == "T":
        maxF = 1.0
        minF = 0.0
        plt.xlabel(r"Disky $\leftarrow$ $T$ (\textonehalf $t_{ff}) \rightarrow$ Filamentary")
    plt.xlim([minF,maxF])
    plt.ylabel("$\%$ TSFE (final)")
    plt.yscale("log")
    plt.savefig("../plots/correlatestructure_"+parameter+"_n"+str(rhofloor).replace(".","p")+"parameter.pdf")

    

if __name__=="__main__":
    run(icsims,"ic")

