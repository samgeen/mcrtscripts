'''
Make images of the cloud/sinks
Sam Geen, June 2016
'''

from startup import *

from pymses.utils import constants as C

import columndensity, sinks, ysos

def FindLscale(boxlen):
    # Find the best match for the length scale bar
    lscales = np.array([1,2,5,10,20,50,100,200,500])
    # Ideal length is around 20% of the box length
    ideal = 0.2*boxlen
    # Find closest value to ideal
    diffs = np.abs(lscales - ideal)
    return int(lscales[diffs == diffs.min()][0])

def ProjectSinks(snap,los):
    up = columndensity.ups[los]
    across = columndensity.acrosses[los]
    sink = sinks.FindSinks(snap)
    sinkx = sink[across]
    sinky = sink[up]
    sinkm = np.atleast_1d(sink.mass)
    return sinkx, sinky, sinkm

# Function "borrowed" from 
# http://www.geophysique.be/2010/11/15/matplotlib-basemap-tutorial-05-adding-some-pie-charts/
def DrawPieChart(ax,X, Y, ratios, size):
    colours = ["red","white"]
    N = len(ratios)
    xy = []
    start = 0.
    first = []
    # Don't draw the line if only one "slice" in the pie
    if max(ratios) != 1.0:
        first = [0]
    for ratio in ratios:
        x = first + np.cos(np.linspace(2*np.pi*start,
                                       2*np.pi*(start+ratio), 30)).tolist()
        y = first + np.sin(np.linspace(2*np.pi*start,
                                       2*np.pi*(start+ratio), 30)).tolist()
        xy1 = zip(x,y)
        xy.append(xy1)
        start += ratio
 
    for i, xyi in enumerate(xy):
        ax.scatter([X],[Y] , marker=(xyi,0), s=size, facecolor=colours[i],
                   lw=0.5)

def MakeImage(snap,los,ax,dolengthscale,ysoage=0.0,timestr=None):
    imtype = "columndensity"
    boxlen = snap.RawData().info["boxlen"]
    numtxt = str(snap.RawData().iout).zfill(5)
    hydro = "NH"
    dmap = columndensity.DensityMap(snap.RawData(),los)
    im = dmap.NH()
    # Plot image map
    yscale = hydrofuncs.yscale(hydro)
    if yscale == "log":
        im = np.log10(im)
    vmin, vmax = hydrofuncs.hydro_range(hydro)
    # Make image mesh
    xl, yl = im.shape
    xarr = np.arange(0,boxlen*1.0000001,boxlen/(xl-1.0))
    yarr = np.arange(0,boxlen*1.0000001,boxlen/(yl-1.0))
    
    # Plot the mesh
    #cmap = "YlGnBu_r"
    cmap = "RdPu_r"
    cax = ax.pcolormesh(xarr,yarr,im,vmin=vmin,vmax=vmax,cmap=cmap)
    cax.set_rasterized(True)
    ax.set_xlim(0,boxlen)
    ax.set_ylim(0,boxlen)

    # Plot sinks
    sinkx, sinky, sinkm = ProjectSinks(snap,los)
    try:
        lenm = len(np.atleast_1d(sinkm))
    except:
        lenm = 0
    if lenm > 0:
        area = np.pi * sinkm / 10.0
        # This is x/y-inverted thanks to the weirdness of 
        #        numpy/pyplot image axes
        if ysoage == 0.0:
            ax.scatter(sinky,sinkx,s=area,c="w",alpha=0.5,
                       edgecolors='w')
        else:
            ysomass = ysos.MassInSnap(snap,los,[ysoage],Aklim=0.0,
                                      persink=True).values()[0]
            frac = ysomass / sinkm
            for x, y, f, a in zip(sinky,sinkx,frac,area):
                r = [f,1-f]
                DrawPieChart(ax,x,y,r,a)

    # Add scale axis
    scalecol = "w"
    if dolengthscale:
        # length scale in pc (hopefully this is what units boxlen is in)
        lscale = FindLscale(boxlen)
        x1 = 0.1 * boxlen
        x2 = x1 + lscale
        y1 = 0.9 * boxlen
        y2 = y1
        ax.plot([x1,x2],[y1,y2],scalecol)
        ax.text(x2,y2, "  "+str(lscale)+" pc",color=scalecol,
                verticalalignment="center",fontsize="x-small")
    # Add times
    if timestr is not None:
        xt = 0.9 * boxlen
        yt = 0.1 * boxlen
        ax.text(xt,yt,timestr,
                horizontalalignment="right",
                verticalalignment="center",
                color=scalecol,fontsize="x-small")
    # Finish up
    ax.set_axis_off()
    ax.set_aspect("equal","datalim")
    ax.axis("off")
    return cax
    fig.subplots_adjust(hspace=0.005, wspace=0.005, 
                        left=0.05,right=0.9,
                        bottom=0.05,top=0.98)
    
    dpi = 200.0
    finches = 4096.0/dpi
    fig.set_size_inches(finches,finches)

def MakeFigure(simname,times,name,los,ysoage=0.0,nonamelabel=False,shape=None):
    times = np.array(times)+0.0 # Paranoid about messing with the inputs
    hydro = "NH"
    ntimes = len(times)
    if shape is not None:
        nrows, ncols = shape
    else:
        nrows = len(times)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, 
                             sharex=True, sharey=True,
                             frameon=False)
    dpi = 200.0
    finches = 512.0/dpi
    # Run for all sims
    dolengthscale = True
    sim = Hamu.Simulation(simname)
    utime = sim.Snapshots()[0].RawData().info["unit_time"].express(C.Myr)
    doxHII = not "-NRT" in simname
    for itime in range(0,ntimes):
        timestr = None
        axtime = times[itime]
        timestr = "{0:.1f}".format(axtime)+" Myr"
        time = times[itime]/utime
        ysos.sim = sim
        snap = sim.FindAtTime(time)
        if ntimes > 1:
            ax = axes.flatten()[itime]
        else:
            ax = axes
        im = MakeImage(snap,los,ax,dolengthscale,ysoage=ysoage,timestr=timestr)
        dolengthscale = False
    # Add colour bar for hydro variable
    print "Making colour bar..."
    cax = fig.add_axes([0.03, 0.1, 0.03, 0.8])
    cbar = fig.colorbar(im,cax=cax)
    label = hydrofuncs.hydro_label((hydro))
    if not "xH" in hydro:
        label = "log("+label+")"
    cbar.set_label(label,fontsize="xx-small",color="k")
    cbar.ax.tick_params(labelsize="x-small",labelcolor="k")
    cbar.solids.set_edgecolor("face")
    # Label limits
    sxl = 0.05
    sxh = 0.9
    sxr = (sxh - sxl) / float(nrows)
    tyl = 0.05
    tyh = 0.98
    tyr = (tyh - tyl) / float(ntimes)
    # Add text - Times
    #for itime2 in range(0,ntimes):
    #    itime = ntimes - itime2 - 1
    #    tax = fig.add_axes([0.01,tyl+tyr*itime2, 0.03, tyr])
    #    axtime = times[itime]
    #    tprefix = ""
    #    timestr = tprefix+"{0:.1f}".format(axtime)+" Myr"
    #    tax.text(0.5,0.5,
    #             timestr,rotation="vertical",fontsize="x-small",
    #             horizontalalignment='center',
    #             verticalalignment='center')
    #    tax.set_axis_off()
    # Add text - Simnames
    if not nonamelabel:
        sax = fig.add_axes([sxl+sxr*isim, 0.01, sxr, 0.03])
        simstr = sim.Label()
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
    #wsname = Hamu.SimData.Settings.Settings()["CurrentWorkspace"]
    folder = "../plots/vis/multiray/"
    MakeDirs(folder)
    figname = folder+"multiray_"+hydro+suffix+".pdf"
    print "Saving figure "+figname+"..."
    fig.subplots_adjust(hspace=0.0, wspace=0.0, 
                        left=0.20,right=1.0,
                        bottom=0.00,top=1.0)
    fig.set_size_inches(finches*nrows/0.8,finches*ncols)
    fig.savefig(figname,
                pad_inches=0,
                dpi=dpi)
    # Crop out borders
    os.system("pdfcrop "+figname+" "+figname)
    print "Done!"

if __name__=="__main__":
    Hamu.Workspace("HIISFE")
    # For Antoine
    simname = "M-RT"
    times = [0.189844823423750E-01,0.271557637924910E-01,
             0.342445791486601E-01,0.414789076094076E-01]
    MakeFigure(simname,times,name="forantoine",los='z')
    '''
    simnames = ["L-RT"]
    tff = 0.05242385002
    times = np.arange(1.0,3.00001,1.0)*tff
    MakeFigure(simnames,times,name="L-RT",los='z')
    simnames = ["M-RT"]
    tff = 0.0155329926
    times = np.arange(1.0,3.00001,1.0)*tff
    MakeFigure(simnames,times,name="M-RT",los='z')
    '''
