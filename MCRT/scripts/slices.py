#import yt
#yt.frontends.ramses.RAMSESDataset._skip_cache=True
import matplotlib as mpl
mpl.use('Agg')

import Hamu

from multiprocessing import Process
import numpy as np
import os, getopt,sys,gc, time
import matplotlib.pyplot as plt

zoom = False



def pymses_func(ro, hydro):
    from pymses.utils import constants as C
    import pymses.analysis.visualization as v
    scop = v.ScalarOperator
    if hydro == "rho":
        unit = ro.info["unit_density"].express(C.H_cc)
        return scop(lambda dset: dset["rho"]*unit)
    if hydro == "P":
        unit = ro.info["unit_pressure"].express(C.barye)
        return scop(lambda dset: dset["P"]*unit)
    if hydro == "T":
        unit = ro.info["unit_temperature"].express(C.K)
        return scop(lambda dset: dset["P"]/dset["rho"]*unit)
    if "xH" in hydro:
        unit = 1.0
        return scop(lambda dset: dset[hydro]*unit)
    if hydro == "Bmag":
        def bmagfunc(dset):
            b = 0.5*(dset["B-left"]+dset["B-right"])
            # Magnitude of the 3-vector for each cell
            return np.sqrt((b**2).sum(axis=1))
        return scop(bmagfunc)
    # None of those? Return unitless
    sco = scop(lambda dset: dset[hydro])
    return sco

def hydro_range(hydro):
    if hydro == "rho":
        return (-1,4)
    if hydro == "P":
        return (None, None) # No idea
    if hydro == "T":
        return (0,5)
    if "xH" in hydro:
        return (-6,0)
    if hydro == "gpe":
        return (None, None)
    if hydro == "Bmag":
        return (None, None)

def hydro_label(hydro):
    if hydro == "rho":
        return "Density / atoms/cm$^{3}$"
    if hydro == "P":
        return "Pressure / dyne"
    if hydro == "T":
        return "Temperature / K"
    if "xH" in hydro:
        return "Ionisation Fraction "+hydro
    if hydro == "gpe":
        return "Gravitational Potential Energy"
    if hydro == "Bmag":
        return "|B| (code units)"

# TODO: MAKE SINGLE FUNCTION FOR SINGLE HYDRO, SNAP
# TODO: PLOT TOGETHER WITH MULTIPLE OUTPUTS FOR COMPARISON

def plotforsnap(snap, hydro,simname=None):
    import pymses
    from pymses.sources.ramses.output import *
    from pymses.analysis.visualization import *
    ro = snap
    amr = ro.amr_source(["rho","P","xHII","B-left","B-right"])
    boxlen = ro.info["boxlen"]
    if not zoom:
        size = np.zeros(2)+1.0
        centre = np.zeros(3)+0.5
        cam  = Camera(center=centre, line_of_sight_axis='z', region_size=size, up_vector='y', 
                      map_max_size=1024, log_sensitive=True)
        zoomtxt = ""
    else:
        size = np.zeros(2)+0.005
        centre = np.zeros(3)+0.5
        centre[1] = 0.5 + 0.010
        cam  = Camera(center=centre, line_of_sight_axis='z', region_size=size, up_vector='y', 
                      map_max_size=1024, log_sensitive=True)
        zoomtxt = "zoom_"
        
    
    # Plot hydro variables
    print hydro
    hydro_op = pymses_func(ro,hydro)
    slc = pymses.analysis.visualization.SliceMap(amr,cam, hydro_op, z=0.0)
    print "Made slice"
    # Only plot if simname is a string? Sure, why not.
    if simname is not None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        br = 0.5 * boxlen
        extent=(-br,br,-br,br)
        imrange = hydro_range(hydro)
        d,u = imrange
        #cax = ax.imshow(slc, interpolation='nearest',\
        #                    extent=extent,vmin=d,vmax=u)
        cax = ax.pcolormesh(np.log10(slc),vmin=d,vmax=u)
        # Add colour bar
        cbar = fig.colorbar(cax)
        label = hydro_label((hydro))
        label = "log("+label+")"
        #cbar.set_label(label)
        # Save figure
        numtxt = str(snap.iout).zfill(5)
        figname = "../plots/slices/"+simname+"/"+hydro+\
            "/slice_"+hydro+"_"+numtxt+".png"
        fig.savefig(figname,\
                        bbox_inches='tight', \
                        transparent=True,\
                        pad_inches=0,
                    dpi=300)
    return slc

plotforsnapHamu = Hamu.Algorithm(plotforsnap)

def runprocs(sim, hydros):
    parallel = True
    last = False
    if last:
        parallel = False
    if parallel:
        maxrun = 10
        running = list()
        procs = list()
        done = 0
        for snap in sim.Snapshots():
            for hydro in hydros:
                proc = Process(target=plotforsnapHamu, args=(snap,hydro,sim.Name()))
                procs.append(proc)
        nprocs = len(procs)
        while done < nprocs:
            for proc in procs:
                # Start process if we have space in the queue
                if not proc.is_alive():
                    if not proc in running:
                        if len(running) < maxrun:
                            print "STARTING", len(running)
                            proc.start()
                            running.append(proc)
                    # Process done
                    else:
                        print "TERMINATING"
                        done += 1
                        running.remove(proc)
                        procs.remove(proc)
            time.sleep(1.0)
    else:
        if last:
            snap = sim.Snapshots()[-1]
        for snap in sim.Snapshots():
            for hydro in hydros:
                plotforsnapHamu(snap,hydro,sim.Name())


def Run(sim,lastonly=False):
    os.system("mkdir ../plots")
    os.system("mkdir ../plots/slices")
    os.system("mkdir ../plots/slices/"+sim.Name())
    hydros = ["rho","P","T","xHII","Bmag"]
    for hydro in hydros:
        os.system("mkdir ../plots/slices/"+sim.Name()+"/"+hydro)
    runprocs(sim, hydros)


if __name__=="__main__":
    ns = ["49","48","47","00"]
    bs = ["02","00","04"]
    #cols = ["r","g","b","k"]
    #lines = ["-","--"]
    for b in bs:
        for n in ns:
            simname = "MCRTN"+n+"B"+b
            if n == "49":
                simname += "T3"
            if n == "48":
                simname += "T6"
            sim = Hamu.Simulation(simname)
            Run(sim)
    # Time-limited emission runs
    #stubs = ["49B00T3","48B00T6","49B02T3","48B02T6"]
    #sims = ["MCRTN"+x for x in stubs]
    #for sim in sims:
    #    s = Hamu.Simulation(sim)
    #    Run(s)
