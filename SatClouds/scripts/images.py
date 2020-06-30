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
    lscales = np.array([1,2,5,10,20,50,100,200,500,1e3,5e3,1e4,5e4,1e5,5e5,1e6])
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

def MakeImage(snap,los,ax,dolengthscale,cmap,label=None,dpi=200.0,simname=None,
              centre=[0.5,0.5,0.5],zoom=1.0,
              labelpos=(0.9,0.1),drawtracks=True):
    imtype = "columndensity"
    cax = None
    if snap is not None:
        boxlen = snap.RawData().info["boxlen"]
        boxlen = 12.5e6
        numtxt = str(snap.RawData().iout).zfill(5)
        hydro = "NH"
        dmap = columndensity.DensityMap(snap.RawData(),los,centre=centre,zoom=zoom)
        im = dmap.NH()
        # Plot image map
        yscale = hydrofuncs.yscale(hydro)
        if yscale == "log":
            im = np.log10(im)
            vmin, vmax = hydrofuncs.hydro_range(hydro)
        # ! HACK!
        vmin = 19
        vmax = 22
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
        nosinks = False
        try:
            sinkx, sinky, sinkm = ProjectSinks(snap,los)
        except:
            nosinks=True
        if not nosinks:
            sinkx -= 0.5*(1.0-zoom)*boxlen
            sinky -= 0.5*(1.0-zoom)*boxlen
            try:
                lenm = len(np.atleast_1d(sinkm))
            except:
                lenm = 0
            if lenm > 0:
                area = np.pi * sinkm / 100.0
                # This is x/y-inverted thanks to the weirdness of 
                #        numpy/pyplot image axes
                #if ysoage == 0.0:
                ax.scatter(sinky,sinkx,s=area,c="w",alpha=0.5,edgecolors='c')

        # Draw star tracks
        if simname is not None and drawtracks:
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

def MakeFigure(sims,snapsdict,name,los,nonamelabel=False,shape=None,dpi=200.0,centre=[0.5,0.5,0.5],zoom=1.0,
               timelabels=False,drawtracks=True):
    hydro = "NH"
    nrows = 1
    ncols = 1
    simnames = [sim.Name() for sim in sims]
    # Make a square array of images
    while len(simnames)*len(snapsdict[simnames[0]]) > nrows*ncols:
        nrows += 1
        ncols += 1
    # Do we fill the whole square?
    square = True
    if len(simnames)*len(snapsdict[simnames[0]]) < nrows*ncols:
        square = False
    # Make figure
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, 
                             sharex=True, sharey=True,
                             frameon=False)
    finches = IMSIZE/dpi
    # Run for all sims
    dolengthscale = True
    isim = -1
    try:
        flaxes = axes.flatten()
    except:
        flaxes = [axes]
    for ax in flaxes:
        ax.set_axis_off()
        ax.set_aspect("equal","datalim")
        ax.axis("off")
    labelpos = (0.9,0.1)
    if square:
        labelpos = (0.9,0.9)
    for sim in sims:
        simname = sim.Name()
        cmap = linestyles.ColourMap(simname)
        label = linestyles.Label(simname)
        for snap in snapsdict[simname]:
            isim += 1
            ax = flaxes[isim]
            if timelabels:
                timestr = "%.1f" % snaptime.Myr(snap)
                label = timestr+" Myr"
            im = MakeImage(snap,los,ax,dolengthscale,cmap,
                           label=label,
                           dpi=dpi, simname=simname,centre=centre,zoom=zoom,labelpos=labelpos,drawtracks=drawtracks)
            dolengthscale = False
    # Make empty frames to remove axis scales, etc
    #for irest in range(isim,nrows*ncols):
    #    square = False
    #    empty = MakeImage(None,los,ax,dolengthscale,cmap,label="",dpi=dpi)
    # Add colour bar for hydro variable
    print "Making colour bar..."
    cax = fig.add_axes([0.28, 0.17, 0.64, 0.03])
    cbar = fig.colorbar(im,cax=cax,orientation="horizontal")
    label = hydrofuncs.hydro_label((hydro))
    if not "xH" in hydro:
        label = "log("+label+")"
    if square:
        cbarcol = "w"
    else:
        cbarcol = "k"
    cbar.set_label(label,color=cbarcol)
    cbar.ax.tick_params(labelcolor=cbarcol)
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
                        bottom=0.0,top=1.0)
    fig.set_size_inches(finches*ncols/0.8,finches*nrows)
    print finches*ncols/0.8,finches*nrows
    fig.savefig(figname,
                pad_inches=0,
                dpi=dpi)
    # Crop out borders
    os.system("pdfcrop "+figname+" "+figname)
    print "Done!"

if __name__=="__main__":
    print "Go!"
    for sim in hamusims.itervalues():
        snaps = [sim.Snapshots()[-1]]
        for halo in ["MW","LMC","SMC"]:
            if halo == "MW":
                centre = [0.74372745, 0.31386477, 0.33639288]
                zoom = 0.1 * 0.186 / 12.5 # 500 kpc in boxlen units (boxlen = 12.5 Mpc)
            if halo == "LMC":
                centre = [0.7384616,  0.2984448,  0.33035183]
                zoom = 0.1 * 0.033266526 / 12.5
            if halo == "SMC":
                centre = [0.73871505, 0.29955167, 0.32425034]
                zoom = 0.1 * 0.031961612 /12.5
            for los in "xyz":
                for zoomout in [True, False]:
                    for snap in snaps:
                        i = snap.OutputNumber()-1
                        print "Snap", i+1
                        simsnapdict = {sim.Name():[snap]}
                        zoomtxt = ""
                        zfact = 1.0
                        ztxt = "RGAL"
                        if zoomout:
                            ztxt = "RVIR"
                            zfact = 10.0
                        MakeFigure([sim],simsnapdict,
                                   name="columns"+sim.Name()+halo+los+ztxt+str(i+1).zfill(5),los=los,
                                   centre=centre,zoom=zoom*zfact,
                                   timelabels=True,drawtracks=False)
