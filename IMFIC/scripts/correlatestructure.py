B'''
Correlate the structure of the cloud to the SFE
Sam Geen, January 2018
'''

from startup import *

import pymses
from pymses.filters import CellsToPoints 
from pymses.utils import constants as C

import triaxfinder, tsfe, regression

#parameter = "F" # Form factor a*b/c2
# Options: F, T, L, M, S, C
parameter = "C" # (1-b^2) / (1-g^2)
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
    posns = cells.points[mask,:]*snap.info["unit_length"].express(C.pc)
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
    if parameter == "C":
        # Basic compactness, radial size of dense parts
        return np.sqrt(np.sum(axes**2))
    if parameter == "L":
        # Longest axis
        return axes[2]
    if parameter == "M":
        # Middlest axis
        return axes[1]
    if parameter == "S":
        # Sweet baby axis
        return axes[0]


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
        label = linestyles.Label(simname)
        plt.scatter(F,sfe*100,s=80,marker="o",c=col,edgecolors='k',zorder=2)
        plt.gca().annotate(r".\ "+label, (F,sfe*100),fontsize="x-small",zorder=1)
    tsfes = np.array(tsfes)
    Fs = np.array(Fs)
    if parameter == "F":
        maxF = np.max([np.max(Fs),10.0])
        minF = np.min([np.min(Fs),0.1])
        plt.xlim([minF,maxF])
        plt.xscale("log")
        plt.xlabel(r"Filamentary $\leftarrow$ $F$ (\textonehalf $t_{ff}) \rightarrow$ Disky")
    if parameter == "T":
        maxF = 1.0
        minF = 0.0
        plt.xlim([minF,maxF])
        plt.xlabel(r"Disky $\leftarrow$ $T$ (\textonehalf $t_{ff}) \rightarrow$ Filamentary")
    if parameter == "C":
        plt.xlabel("Mean radius of cloud / pc")
    if parameter == "L":
        plt.xlabel("Length of longest axis / pc")
    if parameter == "M":
        plt.xlabel("Length of middlest axis / pc")
    if parameter == "S":
        plt.xlabel("Length of shortest axis / pc")
    plt.ylabel("$\%$ TSFE (final)")
    plt.yscale("log")
    plotname = "../plots/correlatestructure_"+parameter+"_n"+str(rhofloor).replace(".","p")+"parameter.pdf"
    plt.savefig(plotname)
    regression.writecoeff(np.log10(Fs),np.log10(tsfes*100),plotname.replace(".pdf","_rxy.txt"))
    

if __name__=="__main__":
    run(icsims,"ic")

