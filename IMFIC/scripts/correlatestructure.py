'''
Correlate the structure of the cloud to the SFE
Sam Geen, January 2018
'''

from startup import *

import pymses
from pymses.filters import CellsToPoints 
from pymses.utils import constants as C

import triaxfinder, regression, plotproperties

#parameter = "F" # Form factor a*b/c2
# Options: F, T, L, M, S, C, V
parameter = "V" # (1-b^2) / (1-g^2)
rhofloor = 10.0

tff_fact = None

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
    if parameter == "V":
        # Volume 4.0/3.0*pi*a*b*c
        return 4.0/3.0*np.pi*a*b*c
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
    # Use different time?
    if tff_fact is not None:
        tcode = tffcloud_code * tff_fact
    snap = sim.FindAtTime(tcode)
    a,b,c = findtriax(snap,cutrho=rhofloor)
    F = FormParameter(a,b,c)
    print F
    # Find last TSFE
    snap = sim.Snapshots()[-1]
    sfe = plotproperties.tsfeinsnap(snap)
    col = linestyles.Colour(simname)
    return sfe,F,col

def makeplot(ax,simnames,plotname,yon):
    print "Running for parameter", parameter
    tsfes = []
    Fs = []
    cols = []
    for simname in simnames:
        sfe, F, col = runforsim(simname)
        tsfes.append(sfe)
        Fs.append(F)
        cols.append(col)
        label = linestyles.Label(simname)
        ax.scatter(F,sfe*100,s=80,marker="o",c=col,edgecolors='k',zorder=2)
        ax.annotate(r".\ "+label, (F,sfe*100),fontsize="x-small",zorder=1)
    tsfes = np.array(tsfes)
    Fs = np.array(Fs)
    if parameter == "F":
        maxF = np.max([np.max(Fs),10.0])
        minF = np.min([np.min(Fs),0.1])
        ax.set_xlim([minF,maxF])
        ax.set_xscale("log")
        ax.set_xlabel(r"Filamentary $\leftarrow$ $F$ (\textonehalf $t_{ff}) \rightarrow$ Disky")
    if parameter == "T":
        maxF = 1.0
        minF = 0.0
        ax.set_xlim([minF,maxF])
        ax.set_xlabel(r"Disky $\leftarrow$ Triaxiality $T \rightarrow$ Filamentary")
    if parameter == "C":
        ax.set_xlabel("Mean radius of cloud $R_E$ / pc")
    if parameter == "L":
        ax.set_xlabel("Length of longest axis $L_E$ / pc")
    if parameter == "M":
        ax.set_xlabel("Length of middlest axis $M_E$ / pc")
    if parameter == "S":
        ax.set_xlabel("Length of shortest axis $S_E$ / pc")
    if parameter == "V":
        ax.set_xlabel("Volume of ellipsoid $V_E$ / pc$^3$")
    #ax.set_yscale("log")
    tfftxt = ""
    if tff_fact is not None:
        tfftxt = "_"+str(tff_fact).replace(".","p")+"tff_"
    #plotname = "../plots/correlatestructure_"+parameter+"_n"+str(rhofloor).replace(".","p")+tfftxt+"parameter.pdf"
    #plt.savefig(plotname)
    regression.writecoeff(np.log10(Fs),np.log10(tsfes*100),plotname.replace(".pdf","_rxy.txt"))
    # Set axis junk
    xlim = ax.get_xlim()
    ax.set_xlim([xlim[0],xlim[1]*1.1])
    ticks = range(5,23)
    ticklabels = {tick:" " for tick in ticks}
    ticklabelsempty = {tick:" " for tick in ticks}
    ticklabels[5] = "5\%"
    ticklabels[10] = "10\%"
    ticklabels[15] = "15\%"
    ticklabels[20] = "20\%"
    ax.yaxis.set_major_locator(plt.NullLocator()) 
    ax.yaxis.set_minor_locator(plt.NullLocator())
    if yon:
        ax.set_ylabel("SFE (final)")
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticklabels.values())
    else:
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticklabelsempty.values())

def run(simnames,plotname,params,suffix):
    global tff_fact, parameter, rhofloor
    nplots = len(params)
    ncols = 3
    nrows = nplots/ncols
    if nplots < ncols:
        ncols = nplots
        nrows = 1
    if nrows*ncols < nplots:
        nrows += 1
    fig, axes = plt.subplots(nrows,ncols)
    try:
        dum = axes.shape
    except:
        axes = np.array([axes])
    # Make figures
    irow = 0
    icol = 0
    yon = True
    for ax, param in zip(axes.flatten(),params):
        # SET VALUES
        tff_fact = 0.5
        parameter = param
        rhofloor = 10.0
        # Run a single plot
        makeplot(ax,simnames,plotname,yon)
        icol += 1
        yon = False
        if icol == ncols:
            icol = 0
            irow += 1
            yon = True
    # Tidy up
    fig.subplots_adjust(wspace=0)
    fig.set_size_inches(6*ncols,5*nrows) 
    plotname = "../plots/correlatestructure_"+suffix+"_"+plotname+".pdf"
    fig.savefig(plotname,bbox_inches='tight',dpi=150) 

if __name__=="__main__":
    params = ["L","M","S","T","C","V"]
    run(icsims,"ic",params,"all")

