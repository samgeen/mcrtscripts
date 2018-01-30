#import yt
#yt.frontends.ramses.RAMSESDataset._skip_cache=True
import customplot

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

fiesta = brewer2mpl.get_map('YlOrRd', 'Sequential', 9,reverse=True).mpl_colormap
coldhot = plt.get_cmap('Spectral_r')

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

def pymses_func(ro, hydro, imtype="projection",rhoweighting=True):
    weight = "rho"
    if not rhoweighting:
        weight = hydro
    class topfunc(object):
        # NOTE: WORKS BADLY, CHANGE WITH SOMETHING ELSE
        # Top operator for projections
        def __init__(self,func):
            self._func = func

        def __call__(self,dset):
            return dset[weight]*self._func(dset)
    if hydro == "NH":
        unit = ro.info["unit_density"].express(C.g_cc) * \
            ro.info["unit_length"].express(C.cm)
        func = lambda dset: dset["rho"]*unit
    elif hydro == "rho":
        unit = ro.info["unit_density"].express(C.H_cc)
        func = lambda dset: dset["rho"]*unit
    elif hydro == "P":
        unit = ro.info["unit_pressure"].express(C.barye)
        func = lambda dset: dset["P"]*unit
    elif hydro == "T":
        mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                     0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
        unit = ro.info["unit_temperature"].express(C.K)
        func = lambda dset: dset["P"]/(dset["rho"]+1e-6)*unit*mufunc(dset)
    elif "xH" in hydro:
        unit = 1.0
        func = lambda dset: dset[hydro]*unit
    elif hydro == "Bmag":
        def bmagfunc(dset):
            b = 0.5*(dset["B-left"]+dset["B-right"])
            # Magnitude of the 3-vector for each cell
            return np.sqrt((b**2).sum(axis=1))
        func = bmagfunc
    # None of those? Return unitless
    else:
        func = lambda dset: dset[hydro]
    # Do projection? Weight by density
    if imtype == "projection" and weight == "rho":
        op = v.FractionOperator(topfunc(func),
                                lambda dset: dset[weight])
    else:
        op = scop(func)
    return op

def hydro_range(hydro):
    if hydro == "rho":
        return (0,8)
    if hydro == "NH":
        return (20.5,25.0)
    if hydro == "P":
        return (None, None) # No idea
    if hydro == "T":
        return (0,5)
    if "xH" in hydro:
        return (0,1)
    if hydro == "gpe":
        return (None, None)
    if hydro == "Bmag":
        return (None, None)

def hydro_label(hydro):
    if hydro == "rho":
        return "Density / atoms/cm$^{3}$"
    if hydro == "NH":
        return "$\mathrm{N}_{\mathrm{H}}$ / cm$^{-2}$"
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

def plotforsnap(snap, ax, hydro,dolengthscale=False,doxHII=True,imtype = "projection"):
    mag = True
    docolumn = imtype == "column"
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
    mapsize = 512

    allfields = ["rho","P","xHII","B-left","B-right","xHeII","xHeIII"]

    def makeray(snap,hydro,mapsize=512,rhoweighting=True,dosurf=docolumn):
        cam  = Camera(center=centre, line_of_sight_axis='z', 
                      region_size=size, up_vector='y', 
                      map_max_size=mapsize, log_sensitive=True)
        hydro_op = pymses_func(snap,hydro,imtype,rhoweighting=rhoweighting)
        rt = pymses.analysis.visualization.raytracing.RayTracer(snap,allfields)
        slc = rt.process(hydro_op,cam,surf_qty=dosurf)
        print "Made "+imtype+" (min/max:", slc.min(), slc.max(), ")"        
        return slc

    def makeslice(snap,hydro,mapsize=512,rhoweighting=True,dosurf=False):
        cam  = Camera(center=centre, line_of_sight_axis='z', 
                      region_size=size, up_vector='y', 
                      map_max_size=mapsize, log_sensitive=True)
        hydro_op = pymses_func(snap,hydro,"slice")
        amr = snap.amr_source(allfields)
        slc = pymses.analysis.visualization.SliceMap(amr,cam, hydro_op, z=0.0)
        print "Made "+imtype+" (min/max:", slc.min(), slc.max(), ")"        
        return slc

    if hydro == "rho" or hydro == "NH":
        cmap = fiesta
    else:
        cmap = coldhot

    if imtype == "projection":
        makeimageHamu = Hamu.Algorithm(makeray)
    elif imtype == "column":
        makeimageHamu = Hamu.Algorithm(makeray)
    else:
        makeimageHamu = Hamu.Algorithm(makeslice)

    # Plot hydro variables
    hydros = [hydro] # HACK, EH
    for hydro in hydros:
        # Get data
        if hydro != "RGB":
            slc = makeimageHamu(snap,hydro,mapsize=mapsize,dosurf=docolumn)
            if hydro == "NH":
                slc /= (1.66e-24/0.74) # Convert from mass to NH
            if not "xH" in hydro:
                slc = np.log10(slc)
        else:
            # Do RGB image with R=T,g=rho,B=xHII
            r = makeimageHamu(snap,"T",mapsize=mapsize)
            g = makeimageHamu(snap,"rho",mapsize=mapsize)
            b = makeimageHamu(snap,"xHII",mapsize=mapsize)
            #import pdb; pdb.set_trace()
            r = np.log10(r)
            g = np.log10(g)
            r = np.flipud(r)
            g = np.flipud(g)
            b = np.flipud(b)
            r = (r-r.min())/(r.max()-r.min())
            g = (g-g.min())/(g.max()-g.min())
            b = (b-b.min())/(b.max()-b.min())
            slc = np.zeros((r.shape[0],r.shape[1],3))
            slc[:,:,0] = r
            slc[:,:,2] = g # Hack
        # Make figure
        br = 0.5 * boxlen
        extent=(-br,br,-br,br)
        #cax = ax.imshow(slc, interpolation='nearest',\
        #                    extent=extent,vmin=d,vmax=u)
        xl = slc.shape[0]
        yl = slc.shape[1]
        xarr = np.arange(-br,br*1.0000001,br*2/(xl-1.0))
        yarr = np.arange(-br,br*1.0000001,br*2/(yl-1.0))
        if hydro != "RGB":
            imrange = hydro_range(hydro)
            d,u = imrange
            cax = ax.pcolormesh(xarr,yarr,slc,vmin=d,vmax=u,cmap=cmap)
        else:
            cax = ax.imshow(slc,extent=extent)
        cax.set_rasterized(True)
        ax.set_xlim(-br,br)
        ax.set_ylim(-br,br)
        # Add point at location of source
        if doxHII:
            ax.plot([0],[0],'o',color='r',markeredgecolor="w",markersize=3.0)
        # Add xHII contours
        hcol = np.array([0.0,1.0,1.0])
        icol = 0.0
        #utime = sim.Snapshots()[0].RawData().info["unit_time"].express(C.Myr)
        #times = np.arange(0.2,2.00001,0.2)/utime
        xHII = makeimageHamu(snap,"xHII",mapsize=mapsize,rhoweighting=False)
        if doxHII:
            try:
                xHIIlim = 100.0*xHII.min()
                col = hcol
                colt = [(col[0],col[1],col[2])]
                ax.contour(xarr,yarr,xHII,[xHIIlim],colors=colt)
            except:
                pass
        # Add scale axis
        if dolengthscale:
            lscale = 10 # pc (hopefully this is what units boxlen is in)
                       #    (if not, oops)
            x1 = 0.1 * 2.0*br - br
            x2 = x1 + lscale
            y1 = 0.1 * 2.0*br - br
            y2 = y1
            scalecol = "w"
            ax.plot([x1,x2],[y1,y2],scalecol)
            ax.text(x2,y2, " "+str(lscale)+" pc",color=scalecol,
                    verticalalignment="center")
        ax.set_axis_off()
        # Add colour bar - xHII contour times
        #axtimes = fig.add_axes([0.15, 0.82, 0.56, 0.05])
        #cdict1 = {'red': ((0.0, 0.0, 0.0),
        #                  (1.0, 0.0, 0.0)),
        #          
        #          'green': ((0.0, 0.0, 0.0),
        #                    (1.0, 1.0, 1.0)),
        #          
        #          
        #          'blue':  ((0.0, 0.0, 0.0),
        #                    (1.0, 1.0, 1.0))
        #          }
        #blackcyan = LinearSegmentedColormap('BlackCyan', cdict1)
        #norm = mpl.colors.Normalize(vmin=0.0, vmax=2.0)
        #cb1 = mpl.colorbar.ColorbarBase(axtimes, cmap=blackcyan,
        #                                norm=norm,
        #                                orientation='horizontal')
        #cb1.set_label("Time / Myr")
        ax.set_aspect("equal","datalim")
        #ax.autoscale_view('tight')
        ax.axis("off")
        return cax

def PlotForSims(simnames, times,hydro="rho", name = "",imtype="projection",
                tstart=(0.0,"")):
    nsim = len(simnames)
    ntimes = len(times)
    fig, axes = plt.subplots(nrows=ntimes, ncols=nsim, sharex=True, sharey=True)
    dpi = 200.0
    finches = 512.0/dpi
    #fig = plt.figure()
    #axes = fig.add_subplot(111)
    # Run for all sims
    dolengthscale = True
    sims = [Hamu.Simulation(simname) for simname in simnames]
    for isim in range(0,nsim):
        simname = simnames[isim]
        sim = sims[isim]
        utime = sim.Snapshots()[0].RawData().info["unit_time"].express(C.Myr)
        doxHII = not "N00" in simname
        for itime in range(0,ntimes):
            time = times[itime]
            # Scale by freefall time for compact cases
            if "_C2" in simname:
                time *= 0.75**3
            elif "_C" in simname:
                time *= 0.5**3
            snap = sim.FindAtTime(time/utime)
            if ntimes > 1 and nsim > 1:
                ax = axes[itime,isim]
            elif nsim > 1:
                ax = axes[isim]
            elif ntimes > 1:
                ax = axes[itime]
            else:
                ax = axes
            im = plotforsnap(snap,ax,hydro, 
                             dolengthscale=dolengthscale,
                             doxHII=doxHII,imtype=imtype)
            dolengthscale = False
    # Add colour bar for hydro variable
    print "Making colour bar..."
    if hydro != "RGB":
        cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(im,cax=cax)

        label = hydro_label((hydro))
        if not "xH" in hydro:
            label = "log("+label+")"
        cbar.set_label(label,fontsize="x-small")
    # Label limits
    sxl = 0.05
    sxh = 0.9
    sxr = (sxh - sxl) / float(nsim)
    tyl = 0.05
    tyh = 0.98
    tyr = (tyh - tyl) / float(ntimes)
    # Add text - Times
    for itime2 in range(0,ntimes):
        itime = ntimes - itime2 - 1
        tax = fig.add_axes([0.01,tyl+tyr*itime2, 0.03, tyr])
        axtime = times[itime] - tstart[0]
        tprefix = ""
        if len(tstart[1]) > 0:
            tprefix = tstart[1] + " + "
        timestr = tprefix+"{0:.1f}".format(axtime)+" Myr"
        tax.text(0.5,0.5,
                 timestr,rotation="vertical",fontsize="x-small",
                 horizontalalignment='center',
                 verticalalignment='center')
        tax.set_axis_off()
    # Add text - Simnames
    for isim in range(0,nsim):
        sax = fig.add_axes([sxl+sxr*isim, 0.01, sxr, 0.03])
        simstr = sims[isim].Label()
        # HACK - replace labels to save space
        simstr = simstr.replace("SN every 0.1 Myr",r"10$ \times $ SNe")
        sax.text(0.5,0.5,
                 simstr,fontsize="x-small",
                 horizontalalignment='center',
                 verticalalignment='center')
        sax.set_axis_off()
    # Save figure
    suffix = ""
    if len(name) > 0:
        suffix = "_"+name
    wsname = Hamu.SimData.Settings.Settings()["CurrentWorkspace"]
    figname = "../plots/vis/multiray/multiray_"+hydro+suffix+"_"+wsname+".pdf"
    print "Saving figure "+figname+"..."
    fig.subplots_adjust(hspace=0.005, wspace=0.005, 
                        left=0.05,right=0.9,
                        bottom=0.05,top=0.98)
    fig.set_size_inches(finches*nsim,finches*ntimes)
    fig.savefig(figname,
                #transparent=True,
                pad_inches=0,
                dpi=dpi)
    print "Done!"

if __name__=="__main__":

    Hamu.Workspace("HIISN")
    simnames = ["N00-SN","N49-SN","N50-SN","N51-SN"]
    times = [5.52,5.52+1.0]
    PlotForSims(simnames,times,name="emissions")
    '''
    simnames = ["N00-NSN","N50-SN","N50-HN"]
    times = [5.52,5.52+1.0]
    PlotForSims(simnames,times,name="hypernova")
    '''
    '''
    Hamu.Workspace("HIISN4")
    simnames = ["N00-NSN","N47-SN","N48-SN","N49-SN"]
    times = [4.25,5.25]
    PlotForSims(simnames,times,name="emissions")
    '''
