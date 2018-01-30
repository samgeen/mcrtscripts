#import yt
#yt.frontends.ramses.RAMSESDataset._skip_cache=True
import matplotlib as mpl
mpl.use('Agg')

from multiprocessing import Process
import numpy as np
import os, getopt,sys,gc, time
import matplotlib.pyplot as plt

zoom = False

import prettyplotlib as ppl
from prettyplotlib import plt
from prettyplotlib import brewer2mpl
import string
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap

red_purple = brewer2mpl.get_map('RdPu', 'Sequential', 9,reverse=True).mpl_colormap

import Hamu

# OK I need to fix this awfulness
import pymses
from pymses.sources.ramses.output import *
from pymses.analysis.visualization import *
from pymses.utils import constants as C
import pymses.analysis.visualization as v
scop = v.ScalarOperator


def runprocs(snaps, hydros,sim):
    parallel = True
    last = False
    if last:
        parallel = False
    if parallel:
        maxrun = 10
        running = list()
        procs = list()
        done = 0
        for snap in snaps:
            # One for slices, one for ray-traced projections
            proc = Process(target=plotforsnap, args=(snap,hydros,sim))
            procs.append(proc)
            #proc = Process(target=plotforsnap, args=(out,hydros,simname,False))
            #procs.append(proc)
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
            outs = [outs[-1]]
        for snap in snaps:
            plotforsnap(snap,hydros,sim)

def pymses_func(ro, hydro):
    if hydro == "rho":
        unit = ro.info["unit_density"].express(C.H_cc)
        #return scop(lambda dset: dset["rho"]*unit)
        return v.FractionOperator(lambda dset: dset["rho"]**2*unit**2,
                                  lambda dset: dset["rho"]*unit)
    if hydro == "P":
        unit = ro.info["unit_pressure"].express(C.barye)
        return scop(lambda dset: dset["P"]*unit)
    if hydro == "T":
        mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                     0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
        unit = ro.info["unit_temperature"].express(C.K)
        return scop(lambda dset: dset["P"]/dset["rho"]*unit*mufunc(dset))
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
        return (0,8)
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

def plotforsnap(snap, hydros,sim):
    mag = True
    imtype = "projection"
    if not mag:
        pymses.sources.ramses.output.RamsesOutput.amr_field_descrs_by_file = \
            {
            "3D": {"hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2, 3]), Scalar("P", 4),
                               Scalar("xHII",5), Scalar("xHeII",6), Scalar("xHeIII",7)]
                   }
            }
    else:
        pymses.sources.ramses.output.RamsesOutput.amr_field_descrs_by_file = \
            {
            "3D": {"hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2, 3]), 
                               Vector("B-left", [4, 5, 6]), 
                               Vector("B-right", [7, 8, 9]), 
                               Scalar("P", 10),
                               Scalar("xHII",11), Scalar("xHeII",12), Scalar("xHeIII",13)]
                   }
            }
        
    ro = snap.RawData()
    boxlen = ro.info["boxlen"]
    size = np.zeros(2)+1.0
    centre = np.zeros(3)+0.5
    cam  = Camera(center=centre, line_of_sight_axis='z', 
                  region_size=size, up_vector='y', 
                  map_max_size=1024, log_sensitive=True)

    def makeray(snap,hydro):
        hydro_op = pymses_func(snap,hydro)
        rt = pymses.analysis.visualization.raytracing.RayTracer(snap,[hydro])
        slc = rt.process(hydro_op,cam)
        print "Made "+imtype+" (min/max:", slc.min(), slc.max(), ")"        
        return slc
    makerayHamu = Hamu.Algorithm(makeray)

    # Plot hydro variables
    for hydro in hydros:
        # Get data
        slc = makerayHamu(snap,"rho")
        slc = np.log10(slc)
        # Make figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        br = 0.5 * boxlen
        extent=(-br,br,-br,br)
        imrange = hydro_range(hydro)
        d,u = imrange
        #cax = ax.imshow(slc, interpolation='nearest',\
        #                    extent=extent,vmin=d,vmax=u)
        xl, yl = slc.shape
        xarr = np.arange(-br,br*1.0000001,br*2/(xl-1.0))
        yarr = np.arange(-br,br*1.0000001,br*2/(yl-1.0))
        cax = ax.pcolormesh(xarr,yarr,slc,vmin=d,vmax=u,cmap=red_purple)
        # Add xHII contours
        hcol = np.array([0.0,1.0,1.0])
        icol = 0.0
        utime = sim.Snapshots()[0].RawData().info["unit_time"].express(C.Myr)
        times = np.arange(0.2,2.00001,0.2)/utime
        ncol = len(times)
        # HACK - JUST DO ONE
        for time in [snap.Time()]:
            #snap = sim.FindAtTime(time)
            col = hcol # * icol/ncol
            colt = [(col[0],col[1],col[2])]
            print colt
            xHII = makerayHamu(snap,"xHII")
            xHIImin = np.min(xHII)
            xHIImax = np.max(xHII)
            try:
                print "CONTOURING", xHIImin, xHIImax
                xHIIlim = xHIImin*10.0
                #import pdb; pdb.set_trace()
                ax.contour(xarr,yarr,xHII,[xHIIlim],colors=colt)
            except:
                pass
            icol += 1
        # Add scale axis
        lscale = 1 # pc (hopefully this is what units boxlen is in)
                   #    (if not, oops)
        x1 = 0.1 * 2.0*br - br
        x2 = x1 + lscale
        y1 = 0.1 * 2.0*br - br
        y2 = y1
        scalecol = "w"
        ax.plot([x1,x2],[y1,y2],scalecol)
        ax.text(x2,y2, " "+str(lscale)+" pc",color=scalecol,
                verticalalignment="center")
        ax.autoscale(tight=True)
        ax.set_axis_off()
        # Add colour bar - Density
        cbar = fig.colorbar(cax)
        label = hydro_label((hydro))
        label = "log("+label+")"
        cbar.set_label(label)
        # Add colour bar - xHII contour times
        '''
        axtimes = fig.add_axes([0.15, 0.82, 0.56, 0.05])
        cdict1 = {'red': ((0.0, 0.0, 0.0),
                          (1.0, 0.0, 0.0)),
                  
                  'green': ((0.0, 0.0, 0.0),
                            (1.0, 1.0, 1.0)),
                  
                  
                  'blue':  ((0.0, 0.0, 0.0),
                            (1.0, 1.0, 1.0))
                  }
        blackcyan = LinearSegmentedColormap('BlackCyan', cdict1)
        norm = mpl.colors.Normalize(vmin=0.0, vmax=2.0)
        cb1 = mpl.colorbar.ColorbarBase(axtimes, cmap=blackcyan,
                                        norm=norm,
                                        orientation='horizontal')
        cb1.set_label("Time / Myr")
        '''
        # Save figure
        simname = sim.Name()
        outnum = str(snap.OutputNumber()).zfill(5)
        figname = "../plots/vis/"+simname+"/hiicontours/hiicontour_"+\
            outnum+".png"
        fig.savefig(figname,\
               bbox_inches='tight', \
               transparent=True,\
               pad_inches=0,
                    dpi=300)

def Run(simname,lastonly=False,num="",hydros=[]):
    os.system("mkdir ../plots/vis/"+simname)
    root = "../plots/vis/"+simname+"/hiicontours/"
    os.system("mkdir "+root)
    if len(hydros) == 0:
        hydros = ["rho","P","T","Bmag"]
    #hydros.append("xHII")
    #for hydro in hydros:
    #    os.system("mkdir "+root+hydro)
    sim = Hamu.Simulation(simname)
    if not lastonly:
        snaps = sim.Snapshots()
    else:
        snaps = [sim.Snapshots()[-1]]
    for snap in snaps:
        plotforsnap(snap, hydros,sim)

if __name__=="__main__":
    lastonly = False
    Run(simname="N48_M4_B02_C",lastonly=lastonly,hydros=["rho"])#$["T","rho","xHII","P"])
