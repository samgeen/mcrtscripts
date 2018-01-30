'''
Plot slices side by side to compare them
Sam Geen, Jul 2014
'''
import customplot
import matplotlib as mpl
mpl.use('Agg')

import os
import Hamu
Hamu.Workspace("MCRT")
import slices
import numpy as np

import matplotlib.pyplot as plt
from pymses.utils import constants as C

def findsnaps(sims,timeToSample):
    '''
    Find the snapshots from a list of simulations at a given time (in Myr)
    '''
    snaps = list()
    for sim in sims:
        # Get time units
        firstSnap = sim.Snapshots()[0].RawData()
        unit_time = firstSnap.info["unit_time"].express(C.Myr)
        # Get snapshot at time
        snap = sim.FindAtTime(timeToSample/unit_time)
        snaps.append(snap)
    return snaps

def imageattime(sims, time, hydro, name):
    '''
    Make an image at a given time (in Myr)
    sims - sims to plot
    time - time to use
    hydro - hydro variable to plot
    name - name to give image
    '''
    sinMyr = 3.15569e13
    snaps = findsnaps(sims,time)
    i = 0
    j = 0
    ncol = 1
    if len(sims) > 3:
        ncol = 2
    for snap, sim in zip(snaps,sims):
        slc = slices.plotforsnapHamu(snap,hydro)
        #slc = np.log10(slc)
        if i == 0 and j == 0:
            image = np.tile(slc,(ncol,len(sims)/ncol))
            snx, sny = slc.shape
            boxlen = snap.RawData().info["boxlen"]
            boxlen *= snap.RawData().info["unit_length"].express(C.pc)
            #time = snap.Time()*snap.RawData().info["unit_time"].express(C.Myr)
        else:
            image[j*snx:(j+1)*snx,i*sny:(i+1)*sny] = slc
        i += 1
        if i >= ncol:
            i = 0
            j += 1
    plt.clf()
    print "Making figure..."
    fig = plt.figure()
    ax = fig.add_subplot(111)
    br = 0.5 * boxlen
    extent=(-br,br,-br,br)
    imrange = slices.hydro_range(hydro)
    nx, ny = image.shape
    d,u = imrange
        #cax = ax.imshow(slc, interpolation='nearest',\
        #                    extent=extent,vmin=d,vmax=u)
    
    cax = ax.pcolormesh(np.log10(image),vmin=d,vmax=u)
    ax.set_xlim([0,ny]) # FU pyplot
    ax.set_ylim([0,nx]) #  "     "
    i = 0
    j = 0
    for snap, sim in zip(snaps,sims):
        ax.text(10+snx*i,10+sny*j,sim.Label(),color="w")
        i += 1
        if i >= ncol:
            i = 0
            j += 1
    # TIME TIME TIIIIIME
    ax.text(10,sny*2-100,"{0:0.1f}".format(time-1.25)+" Myr",
            color="w",weight="bold")
    # Add colour bar
    cbar = fig.colorbar(cax)
    label = slices.hydro_label((hydro))
    label = "log("+label+")"
    cbar.set_label(label)
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().get_yaxis().set_ticks([])
    ax.set_aspect("equal")
    fig.savefig("../plots/vis/compare/slices/"+name+".png",density=200,bbox_inches="tight")
    plt.close(fig)
    print "Done figure"


if __name__=='__main__':
    
    simnames = ["N00B00","N00B02"]
    ns = ["00","47","48","49"]
    bs = ["02"]
    for b in bs:
        try:
            os.mkdir("../plots/slices/compare/B"+b)
        except:
            pass
        for i in np.arange(0.0,100.0):
            numtxt = str(int(i)).zfill(5)
            time = 10.0 * i / 100.0
            sims = [Hamu.Simulation("N"+n+"_M4_B"+b) for n in ns]
            imageattime(sims,time,"rho","B"+b+"/compare"+numtxt)
        os.system("convert "+\
                  "../plots/vis/compare/slices/B"+b+"/compare*.png "+\
                  "../plots/vis/compare/slices/compare_B"+b+".gif")
