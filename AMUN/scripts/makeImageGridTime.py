'''
Make images of the cloud/sinks
Rebekka Bieri and Sam Geen, August 2018
'''

from startup import *

from pymses.utils import constants as C

import columndensity, rayMap, sliceMap, sinks, ysos

from matplotlib import rc

rc('axes', labelsize=8)
rc('axes', labelsize=8,linewidth=0.5)
rc('xtick', labelsize=8)
rc('ytick', labelsize=8)

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

def _createColDensMap(snap,los=None,hydro='rho',wsink=False,zoom=1.0):
    dmap = columndensity.DensityMap(snap,los,zoom=zoom)
    im = dmap.NH()
    boxlen = snap.info["boxlen"]
    return [im, boxlen*zoom]

def _createColDensMap_sink(snap,los=None,hydro='rho',wsink=True,zoom=1.0):
    dmap = columndensity.DensityMap(snap,los,zoom=zoom)
    im = dmap.NH()
    boxlen = snap.info["boxlen"]
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(snap,los)
    # Shift sink position if zooming in
    sinkx -= 0.5*(1.0-zoom)*boxlen
    sinky -= 0.5*(1.0-zoom)*boxlen
    return [im, sinkx, sinky, sinkm, boxlen*zoom]

def _createRayMap(snap,los=None,hydro='rho',wsink=False,zoom=1.0):
    dmap = rayMap.RayTraceMap(snap,hydro,los,zoom=zoom)
    im   = dmap.getRaytraceMap()
    boxlen = snap.info["boxlen"]
    return [im, boxlen*zoom]

def _createRayMap_sink(snap,los=None,hydro='rho',wsink=True,zoom=1.0):
    ro = snap.RawData()
    dmap = rayMap.RayTraceMap(snap,hydro,los,zoom=zoom)
    im   = dmap.getRaytraceMap()
    boxlen = ro.info["boxlen"]
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(ro,los)
    # Shift sink position if zooming in
    sinkx -= 0.5*(1.0-zoom)*boxlen
    sinky -= 0.5*(1.0-zoom)*boxlen
    return [im, sinkx, sinky, sinkm, boxlen*zoom]

def _createSliceMapStar(snap,hydro='rho',los=None,zoom=1.0):
    dmap = sliceMap.SliceMap(snap,hydro,los,starC=True,zoom=zoom)
    im   = dmap.getSliceMap()
    boxlen = snap.info["boxlen"] 
    return [im, boxlen*zoom]

def _createSliceMap(snap,hydro='rho',los=None,zoom=1.0):
    dmap = sliceMap.SliceMap(snap,hydro,los,zoom=zoom)
    im   = dmap.getSliceMap()
    boxlen = snap.info["boxlen"] 
    return [im, boxlen*zoom]

def _createSliceMap_sink(snap,hydro='rho',los=None,zoom=1.0):
    dmap = sliceMap.SliceMap(snap,hydro,los,zoom=zoom)
    im   = dmap.getSliceMap()
    boxlen = snap.info["boxlen"] 
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(snap,los)
    # Shift sink position if zooming in
    sinkx -= 0.5*(1.0-zoom)*boxlen
    sinky -= 0.5*(1.0-zoom)*boxlen
    return [im, sinkx, sinky, sinkm, boxlen*zoom]
 
def MakeImageHamu(data,hydro,wsink,ax,dolengthscale,cmap,plottime=False,timeL=None,label=None):

    if wsink:
        im, sinkx, sinky, sinkm, boxlen = data
    else:
        im, boxlen = data

    # Specific hack for Lcool (ignore heating, cooling = +ve)
    if hydro == "Lcool" or "xrayemission" in hydro or "ionemission" in hydro:
        #im = -im
        im[im <= 0.0] = im[im > 0.0].min()*0.1
        
    # Plot image map
    yscale = hydrofuncs.yscale(hydro)
    vmin, vmax = hydrofuncs.hydro_range(hydro)
    if yscale == "log":
        if np.min(im) > 0.0:
            im = np.log10(im)
        else:
            # Do log scale where -ve values exist
            yscale = "symlog"
            im2 = np.log10(im)
            im2[im < 0.0] = np.log10(-im[im < 0.0])
    # Make image mesh
    xl, yl = im.shape
    xarr = np.arange(0,boxlen*1.0000001,boxlen/(xl-1.0))
    yarr = np.arange(0,boxlen*1.0000001,boxlen/(yl-1.0))
    
    # Plot the mesh
    #cmap = "YlGnBu_r"
    #cmap = "RdPu_r"
    cax = ax.pcolormesh(xarr,yarr,im,vmin=vmin,vmax=vmax,cmap=cmap)
    cax.set_rasterized(True)
    ax.set_xlim(0,boxlen)
    ax.set_ylim(0,boxlen)
   
    if wsink: 
        # Plot sinks
        try:
            lenm = len(np.atleast_1d(sinkm))
        except:
            lenm = 0
        if lenm > 0:
            area = np.pi * sinkm / 50.0
            ax.scatter(sinky,sinkx,s=area,c="w",alpha=0.5,edgecolors='w')

    # Add scale axis
    scalecol = "w"
    if dolengthscale:
        # length scale in pc (hopefully this is what units boxlen is in)
        lscale = FindLscale(boxlen)
        x1 = 0.7 * boxlen 
        x2 = x1 + lscale
        y1 = 0.9 * boxlen
        y2 = y1
        ax.plot([x1,x2],[y1,y2],scalecol)
        ax.text(x2,y2, "  "+str(lscale)+" pc",color=scalecol,
                verticalalignment="center",fontsize="large")
    # Add label
    if label is not None:
        xt = 0.9 * boxlen
        yt = 0.1 * boxlen
        ax.text(xt,yt,label,
                horizontalalignment="right",
                verticalalignment="center",
                color=scalecol,fontsize="large")
    # Add time 
    if plottime:
        xt = 0.2 * boxlen
        yt = 0.9 * boxlen
        ax.text(xt,yt,timeL,
                horizontalalignment="right",
                verticalalignment="center",
                color=scalecol,fontsize="large")
    return cax

def MakeFigure(simnames,times,name,los=None,hydro="rho",Slice=False,wsink=False,starC=False,
                nonamelabel=False,timeL=None,shape=None,dpi=200.0,zoom=1.0):
    ncols = len(simnames)
    nrows = len(times)
 
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, 
                             sharex=True, sharey=True,
                             frameon=False)
    finches = IMSIZE/dpi
    
    axes = np.atleast_1d(axes)

    if len(axes) <= 1:
        finches *= 1.5 # heuristic
    
    # Originally this was run through Hamu, but we don't really need another level of Hamu for otherwise quick functions
    if wsink:
        createColDensMap_sink = _createColDensMap_sink
        createRayMap_sink     = _createRayMap_sink
        createSliceMap_sink   = _createSliceMap_sink
    else:
        createColDensMap   = _createColDensMap
        createRayMap       = _createRayMap
        createSliceMap     = _createSliceMap
        createSliceMapStar = _createSliceMapStar

    # Run for all sims
    dolengthscale = False
    plottime      = True
    isim = -1
    for ax in axes.flatten():
        ax.set_axis_off()
        ax.set_aspect("equal", "datalim")
        ax.axis("off")

    for ii, time in enumerate(times):
        for simname in simnames:
            isim += 1
            cmap  = linestyles.ColourMap(simname, hydro)
            sim   = Hamu.Simulation(simname)
            snap  = sim.FindAtTime(time)
            ax    = axes.flatten()[isim]
          
            if Slice:
                if wsink:
                    data  = createSliceMap_sink(snap,hydro,los,zoom)      
                else:
                    if starC:
                        data  = createSliceMapStar(snap,hydro,los,zoom)     
                    else:
                        data  = createSliceMap(snap,hydro,los,zoom)     
            else: 
                if hydro == 'NH':
                    if wsink: 
                        data  = createColDensMap_sink(snap,los,hydro,wsink,zoom)      
                    else: 
                        data  = createColDensMap(snap,los,hydro,wsink,zoom)       
                else:
                    if wsink: 
                        data = createRayMap_sink(snap,los,hydro,wsink,zoom)    
                    else: 
                        data = createRayMap(snap,los,hydro,wsink,zoom)   
 
            if (simname == simnames[-1] and ii == 0):
                dolengthscale = True 
            if (simname == simnames[0]) and len(axes) > 1:
                plottime      = True
            label =  linestyles.Label(simname)
            if len(axes) == 1:
                label = None
            im    = MakeImageHamu(data,hydro,wsink,ax,dolengthscale,cmap,
                           plottime, timeL[ii],label = label)
            plottime      = False
            dolengthscale = False

    # Add colour bar for hydro variable
    print "Making colour bar..."
    # ---- This cax works well with nrow = 1, may need adjusting, could be improved 
    # Colorbar at the bottom of the plots
    #cax  = fig.add_axes([0.2, -0.022, 0.4, 0.02])
    # Colorbar at the top of all the plots
    cax  = fig.add_axes([0.4, 0.98, 0.4, 0.02])
    cbar = fig.colorbar(im,cax=cax,orientation="horizontal")
    label = hydrofuncs.hydro_label((hydro))
    if not "xH" in hydro:
        label = "log("+label+")"
    #cbar.set_label(label,fontsize="medium",color="k")
    #cbar.ax.tick_params(labelsize="medium",labelcolor="k")
    #cbar.solids.set_edgecolor("face")
    cbar.set_label(label,fontsize="large",color="w")
    cbar.ax.tick_params(labelsize="large",labelcolor="w")
    cbar.solids.set_edgecolor("face")
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')

    # Save figure
    suffix = ""
    if len(name) > 0:
        suffix = "_"+name
    #wsname = Hamu.SimData.Settings.Settings()["CurrentWorkspace"]
    if Slice:
        folder = "../plots/vis/slice/"
        MakeDirs(folder)
        figname = folder+"sliceTime_"+hydro+suffix+".pdf"
    else:
        folder = "../plots/vis/multiray/"
        MakeDirs(folder)
        figname = folder+"multirayTime_"+hydro+suffix+".pdf"
    print "Saving figure "+figname+"..."
    fig.subplots_adjust(hspace=0.0, wspace=0.0, 
                        left=0.20,right=1.0,
                        bottom=0.00,top=1.0)
    #fig.set_size_inches(finches*ncols/0.8,finches*(nrows+1))
    fig.set_size_inches(finches*ncols/0.8,finches*(nrows))
    fig.savefig(figname,
                pad_inches=0,
                dpi=dpi)
    # Crop out borders
    os.system("pdfcrop "+figname+" "+figname)
    print "Done!"

if __name__=="__main__":

    for mass in [30,60,120][::-1]:
        smass = str(mass)
        simset = ["NOFB","UV_"+smass,"UVWIND_"+smass]
        setname = "windset_"+smass+"Msun"
        #times = np.array([0.5, 0.75, 1.])
        times = np.array([0.9]) # 3.5 Myr = tstarformed + 0.2 Myr 
        zoom = 0.5
        setname = setname+str(times[-1])+"tff_"+"zoom"+str(zoom)+"_"
        timeL = [str(x)+r' t$_{ff}$' for x in times]
        timescode = times * tffcloud_code
        for hydro in ["ionemission","xrayemission","xrayemission2"][::-1]:
            MakeFigure([simset[-1]],[timescode[-1]],name=setname+"windonly",los='z',hydro=hydro,Slice=False,wsink=True,
                       timeL=[timeL[-1]],zoom=zoom)
        for hydro in ["Lcool","T","rho","xHII"]:
            MakeFigure([simset[-1]],[timescode[-1]],name=setname+"windonly",los='z',hydro=hydro,Slice=True,wsink=True,
                       timeL=[timeL[-1]],zoom=zoom)
        MakeFigure([simset[-1]],[timescode[-1]],name=setname+"windonly",los='z',hydro='NH',Slice=False,wsink=True,
                   timeL=[timeL[-1]],zoom=zoom)
        MakeFigure(simset,timescode,name=setname,los='z',hydro='NH',Slice=False,wsink=True,timeL=timeL,zoom=zoom)
        MakeFigure(simset,timescode,name=setname,los='z',hydro='T',Slice=True,wsink=True,timeL=timeL,zoom=zoom)
        
