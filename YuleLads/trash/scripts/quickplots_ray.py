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

# All of this is very bad
from pymses.utils import constants as C
import pymses.analysis.visualization as v
scop = v.ScalarOperator
import pymses
from pymses.sources.ramses.output import *
from pymses.analysis.visualization import *

#import cleaninfo

import simlabels

import sinks

red_purple = brewer2mpl.get_map('RdPu', 'Sequential', 9,reverse=True).mpl_colormap

justimage = True

def runprocs(outs, hydros,simname):
    parallel = True
    last = False
    if len(outs) == 1:
        last = True
    if last:
        parallel = False
    if parallel:
        maxrun = 2
        running = list()
        procs = list()
        done = 0
        for out in outs:
            # One for slices, one for ray-traced projections
            proc = Process(target=plotforsnap, args=(out,hydros,simname,True))
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
        for out in outs:
            plotforsnap(out,hydros,simname)

class MaxTempOperator(Operator):
    def __init__(self, ro):
        self._unit = ro.info["unit_temperature"].express(C.K)
        def Tfunc(dset):
            mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                       0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
            T = dset["P"]/dset["rho"]*self._unit*mufunc(dset)
            return T
        d = {"T": Tfunc}
        Operator.__init__(self, d, is_max_alos=True)

    def operation(self, int_dict):
        mapT = int_dict.values()[0]
        return mapT

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
        return MaxTempOperator(ro)
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

def plotforsnap(out, hydros,simname,ray=False):
    mag = True
    imtype = "projection"
    zoomtxt = ""
    out = out.replace("\n","")
    print "Reading file", out
    oloc = out.find("output_")
    num = int(out[oloc+7:oloc+12])
    numtxt = out[oloc+7:oloc+12]
    if zoom:
        zoomtxt = "zoom_"
    newhydros = []
    for hydro in hydros:
        figname = "../plots/"+simname+"/"+imtype+"s/"+hydro+"/"+\
                  imtype+"_"+hydro+"_"+zoomtxt+numtxt+".png"
        if not os.path.isfile(figname):
            newhydros.append(hydro)
    hydros = newhydros
    if len(hydros) == 0:
        return

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
        
    folder = out[0:oloc]
    ro = RamsesOutput(folder, num)
    allhydros = ["rho","P","xHII","B-left","B-right","xHeII","xHeIII"]
    amr = ro.amr_source(allhydros)
    boxlen = ro.info["boxlen"]
    if not zoom:
        size = np.zeros(2)+1.0
        centre = np.zeros(3)+0.5
        cam  = Camera(center=centre, line_of_sight_axis='z', region_size=size, up_vector='y', 
                      map_max_size=512, log_sensitive=True)
    else:
        size = np.zeros(2)+0.005
        centre = np.zeros(3)+0.5
        centre[1] = 0.5 + 0.010
        cam  = Camera(center=centre, line_of_sight_axis='z', region_size=size, up_vector='y', 
                      map_max_size=512, log_sensitive=True)
        

    # Get sink info
    sink = sinks.FindSinks(ro)
        
    # Plot hydro variables
    for hydro in hydros:
        print hydro
        hydro_op = pymses_func(ro,hydro)
        rt = pymses.analysis.visualization.raytracing.RayTracer(ro,allhydros)
        slc = rt.process(hydro_op,cam)
        print "Made "+imtype+" (min/max:", slc.min(), slc.max(), ")"        
        slc = np.log10(slc)
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
        if not justimage:
            cax = ax.pcolormesh(xarr,yarr,slc,vmin=d,vmax=u,cmap=red_purple)
        else:
            cax = ax.imshow(slc,vmin=d,vmax=u,cmap=red_purple)
        # Add scale axis
        if not justimage:
            lscale = 5 # pc (hopefully this is what units boxlen is in)
                       #    (if not, oops)
            x1 = 0.1 * 2.0*br - br
            x2 = x1 + lscale
            y1 = 0.1 * 2.0*br - br
            y2 = y1
            scalecol = "w"
            ax.plot([x1,x2],[y1,y2],scalecol)
            ax.text(x2,y2, " "+str(lscale)+" pc",color=scalecol,
                    verticalalignment="center")
            # Add colour bar
            cbar = fig.colorbar(cax)
            label = hydro_label((hydro))
            label = "log("+label+")"
            cbar.set_label(label)

        #ax.autoscale(tight=True)
        ax.set_axis_off()

        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent = im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
        forceAspect(ax)
        # Make scatter overlaid
        sx = sink.y / (2*br)
        sy = 1.0 - sink.x / (2*br)
        sz = sink.z / (2*br)
        inside = (sx >= 0)*(sx <= 1)*(sy >= 0)*(sy <= 1)*(sz >= 0)*(sz <= 1)
        sx = sx[inside]
        sy = sy[inside]
        sz = sz[inside]
        if len(sx) > 0:
            ext = ax.get_images()[0].get_extent()
            sx = sx * (ext[1] - ext[0]) + ext[0]
            sy = sy * (ext[3] - ext[2]) + ext[2]
            ax.scatter(sx,sy,color="w")

        # Save figure
        figname = "../plots/"+simname+"/"+imtype+"s/"+hydro+"/"+\
                  imtype+"_"+hydro+"_"+zoomtxt+numtxt+".png"

        fig.savefig(figname,
                    #bbox_inches='tight',
                    transparent=True,
                    #pad_inches=0,
                    dpi=150)
            
def Run(folder=".",simname="",lastonly=False,num="",hydros=[]):
    #cleaninfo.run()
    os.system("mkdir ../plots")
    os.system("mkdir ../plots/"+simname)
    #os.system("mkdir ../plots/"+simname+"/slices")
    os.system("mkdir ../plots/"+simname+"/projections")
    #os.system("mkdir plots/projs")
    if len(hydros) == 0:
        hydros = ["rho","P","T","xHII","Bmag"]
    #hydros = ["gpe"]
    if zoom:
        hydros = ["rho"]
    for hydro in hydros:
        #os.system("mkdir ../plots/"+simname+"/slices/"+hydro)
        os.system("mkdir ../plots/"+simname+"/projections/"+hydro)
        #os.system("mkdir plots/projs/"+hydro)
    p = os.popen("ls "+folder+"/output_?????/info_?????.txt")
    outs = p.readlines()
    if lastonly:
        outs = [outs[-1]]
    if num != "":
        print outs, num
        outs = [x for x in outs if num in x]
    runprocs(outs, hydros,simname)
    #for out in outs:
    #    plotforsnap(out, hydros)
    #for hydro in hydros:
    #    os.system("mv *Slice*"+hydro+"*.png plots/slices/"+hydro+"/")
    #    os.system("mv *Projection*"+hydro+"*.png plots/projs/"+hydro+"/")

if __name__=="__main__":
    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv,"l")
    lastonly = False
    for opt, arg in opts:
        if opt == "-l":
            print "READING LAST OUTPUT ONLY"
            lastonly = True
    allsims = simlabels.allsims
    for simname in allsims:
        folder = "/home/sgeen/MC_RT/runs_anais/48_yulelads/"+simname+"/"
        Run(folder=folder,simname=simname,lastonly=lastonly,hydros=["rho"])#,"rho","xHII","P"])
