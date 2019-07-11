'''
Make images of the cloud/sinks
Sam Geen, June 2016
'''

from startup import *

from pymses.utils import constants as C

import columndensity, sinks, ysos, trackstars

IMSIZE = columndensity.IMSIZE

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

def MakeImage(snap,los,ax,dolengthscale,cmap,label=None,dpi=200.0,simname=None,zoom=1.0,labelpos=(0.9,0.1)):
    imtype = "columndensity"
    cax = None
    if snap is not None:
        boxlen = snap.RawData().info["boxlen"]
        numtxt = str(snap.RawData().iout).zfill(5)
        hydro = "NH"
        dmap = columndensity.DensityMap(snap.RawData(),los,zoom=zoom)
        im = dmap.NH()
        # Plot image map
        yscale = hydrofuncs.yscale(hydro)
        if yscale == "log":
            im = np.log10(im)
            vmin, vmax = hydrofuncs.hydro_range(hydro)
        # Make image mesh
        xl, yl = im.shape
        xarr = np.arange(0,zoom*boxlen*1.0000001,zoom*boxlen/(xl-1.0))
        yarr = np.arange(0,zoom*boxlen*1.0000001,zoom*boxlen/(yl-1.0))
    
        # Plot the mesh
        #cmap = "YlGnBu_r"
        #cmap = "RdPu_r"
        cax = ax.pcolormesh(xarr,yarr,im,vmin=vmin,vmax=vmax,cmap=cmap)
        cax.set_rasterized(True)
        ax.set_xlim(0,zoom*boxlen)
        ax.set_ylim(0,zoom*boxlen)

        # Plot sinks
        sinkx, sinky, sinkm = ProjectSinks(snap,los)
        sinkx -= 0.5*(1.0-zoom)*boxlen
        sinky -= 0.5*(1.0-zoom)*boxlen
        try:
            lenm = len(np.atleast_1d(sinkm))
        except:
            lenm = 0
        if lenm > 0:
            area = np.pi * sinkm / 10.0
        # This is x/y-inverted thanks to the weirdness of 
        #        numpy/pyplot image axes
        #if ysoage == 0.0:
        ax.scatter(sinky,sinkx,s=area,c="w",alpha=0.5,edgecolors='w')

        # Draw star tracks
        if simname is not None:
            tracks, tracktimes, dummy = trackstars.runforsim(simname)
            if los == "z":
                xtrack = 0
                ytrack = 1
            if los == "x":
                xtrack = 1
                ytrack = 2
            if los == "y":
                xtrack = 2
                ytrack = 0
            for track in tracks.itervalues():
                track -= 0.5*(1.0-zoom)*boxlen
                # Plot track
                ax.plot(track[:,ytrack],track[:,xtrack],alpha=0.5,color="w")
                x0 = track[-1,ytrack]
                y0 = track[-1,xtrack]
                x1 = x0-track[-2,ytrack]
                y1 = y0-track[-2,xtrack]
                x1 *= 0.01
                y1 *= 0.01
                # Plot arrow at end of track
                ax.arrow(x0, y0, x1, y1, shape='full', 
                         lw=0, length_includes_head=False, head_width=0.02*boxlen*zoom,
                         alpha=0.5,color="w")

        # Add scale axis
        scalecol = "w"
        if dolengthscale:
            # length scale in pc (hopefully this is what units boxlen is in)
            lscale = FindLscale(boxlen*zoom)
            x1 = 0.1 * boxlen*zoom
            x2 = x1 + lscale
            y1 = 0.9 * boxlen*zoom
            y2 = y1
            ax.plot([x1,x2],[y1,y2],scalecol)
            ax.text(x2,y2, "  "+str(lscale)+" pc",color=scalecol,
                    verticalalignment="center",fontsize="x-large")
        # Add label
        if label is not None:
            xt = labelpos[0] * boxlen*zoom
            yt = labelpos[1] * boxlen*zoom
            ax.text(xt,yt,label,
                    horizontalalignment="right",
                    verticalalignment="center",
                    color=scalecol,fontsize="x-large")
    # Finish up
    #ax.set_axis_off()
    #ax.set_aspect("equal","datalim")
    #ax.axis("off")
    return cax

def MakeFigure(simnames,times,name,los,nonamelabel=False,shape=None,dpi=200.0,zoom=1.0,timelabels=False):
    hydro = "NH"
    nrows = 1
    ncols = 1
    # Make a square array of images
    while len(simnames)*len(times) > nrows*ncols:
        nrows += 1
        ncols += 1
    # Do we fill the whole square?
    square = True
    if len(simnames)*len(times) < nrows*ncols:
        square = False
    # Make figure
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, 
                             sharex=True, sharey=True,
                             frameon=False)
    finches = IMSIZE/dpi
    # Run for all sims
    dolengthscale = True
    isim = -1
    for ax in axes.flatten():
        ax.set_axis_off()
        ax.set_aspect("equal","datalim")
        ax.axis("off")
    labelpos = (0.9,0.1)
    if square:
        labelpos = (0.9,0.9)
    for simname in simnames:
        cmap = linestyles.ColourMap(simname)
        sim = Hamu.Simulation(simname)
        label = linestyles.Label(simname)
        for time in times:
            isim += 1
            snap = sim.FindAtTime(time)
            ax = axes.flatten()[isim]
            if timelabels:
                timestr = "%.1f" % snaptime.Myr(snap)
                label = timestr+" Myr"
            im = MakeImage(snap,los,ax,dolengthscale,cmap,
                           label=label,
                           dpi=dpi, simname=simname,zoom=zoom,labelpos=labelpos)
            dolengthscale = False
    # Make empty frames to remove axis scales, etc
    #for irest in range(isim,nrows*ncols):
    #    square = False
    #    empty = MakeImage(None,los,ax,dolengthscale,cmap,label="",dpi=dpi)
    # Add colour bar for hydro variable
    print "Making colour bar..."
    cax = fig.add_axes([0.43, 0.10, 0.54, 0.03])
    cbar = fig.colorbar(im,cax=cax,orientation="horizontal")
    label = hydrofuncs.hydro_label((hydro))
    if not "xH" in hydro:
        label = "log("+label+")"
    if square:
        cbarcol = "w"
    else:
        cbarcol = "k"
    cbar.set_label(label,fontsize="x-large",color=cbarcol)
    cbar.ax.tick_params(labelsize="x-large",labelcolor=cbarcol)
    cbar.solids.set_edgecolor("face")
    # Save figure
    suffix = ""
    if len(name) > 0:
        suffix = "_"+name
    #wsname = Hamu.SimData.Settings.Settings()["CurrentWorkspace"]
    folder = "../plots/vis/multiray/"
    MakeDirs(folder)
    figname = folder+"multiray_"+hydro+suffix+".pdf"
    print "Saving figure "+figname+"..."
    fig.subplots_adjust(hspace=0.02, wspace=0.02, 
                        left=0.20,right=1.0,
                        bottom=0.00,top=1.0)
    fig.set_size_inches(finches*ncols/0.8,finches*nrows)
    fig.savefig(figname,
                pad_inches=0,
                dpi=dpi)
    # Crop out borders
    os.system("pdfcrop "+figname+" "+figname)
    print "Done!"

if __name__=="__main__":
    time = 1.5*tffcloud_code
    times = np.array([1.0,1.5,2.0,2.5])*tffcloud_code
    zoom = 0.7 # 1.0
    for i in range(0,13):
        for los in "xyz":
            MakeFigure([icsims[i]],times,name="drifttimes_ic"+str(i)+los,los=los,zoom=zoom,timelabels=True)
    for los in "xyz":
        MakeFigure([imfsims[2]],times,name="drifttimes"+los,los=los,zoom=zoom,timelabels=True)
    MakeFigure(imfsims,[time],name="imf",los='z',zoom=zoom)
    MakeFigure(icsims,[time],name="ic",los='z',zoom=zoom)
