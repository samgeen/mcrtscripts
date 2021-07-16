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

def _createColDensMap(snap,los=None,hydro='rho',wsink=False):
    dmap = columndensity.DensityMap(snap,los,hydro)
    im = dmap.NH()
    boxlen = snap.info["boxlen"]
    return [im, boxlen]

def _createColDensMap_sink(snap,los=None,hydro='rho',wsink=True):
    dmap = columndensity.DensityMap(snap,los,hydro)
    im = dmap.NH()
    boxlen = snap.info["boxlen"]
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(snap,los)
    return [im, sinkx, sinky, sinkm, boxlen]

def _createRayMap(snap,los=None,hydro='rho',wsink=False):
    dmap = rayMap.RayTraceMap(snap,hydro,los)
    im   = dmap.getRaytraceMap()
    boxlen = snap.info["boxlen"]
    return [im, boxlen]

def _createRayMap_sink(snap,los=None,hydro='rho',wsink=True):
    dmap = rayMap.RayTraceMap(snap,hydro,los)
    im   = dmap.getRaytraceMap()
    boxlen = snap.info["boxlen"]
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(snap,los)
    return [im, sinkx, sinky, sinkm, boxlen]

def _createSliceMapStar(snap,hydro='rho',los=None):
    dmap = sliceMap.SliceMap(snap,hydro,los,starC=True)
    im   = dmap.getSliceMap()
    boxlen = snap.info["boxlen"] 
    return [im, boxlen]

def _createSliceMap(snap,hydro='rho',los=None):
    dmap = sliceMap.SliceMap(snap,hydro,los)
    im   = dmap.getSliceMap()
    boxlen = snap.info["boxlen"] 
    return [im, boxlen]

def _createSliceMap_sink(snap,hydro='rho',los=None):
    dmap = sliceMap.SliceMap(snap,hydro,los)
    im   = dmap.getSliceMap()
    boxlen = snap.info["boxlen"] 
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(snap,los)
    return [im, sinkx, sinky, sinkm, boxlen]
 
def MakeImageHamu(data,hydro,wsink,ax,dolengthscale,cmap,label=None):

    if wsink:
        im, sinkx, sinky, sinkm, boxlen = data
    else:
        im, boxlen = data
   
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
        x1 = 0.1 * boxlen
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
                color=scalecol,fontsize="medium")
    return cax

def MakeFigure(simnames,time,name,los=None,hydro="rho",Slice=False,wsink=False,starC=False,nonamelabel=False,shape=None,dpi=200.0):
    ncols = 4
    nrows = max(int(len(simnames)/ncols),1)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, 
                             sharex=True, sharey=True,
                             frameon=False)
    finches = IMSIZE/dpi

    if wsink:
        createColDensMap_sink = Hamu.Algorithm(_createColDensMap_sink)
        createRayMap_sink     = Hamu.Algorithm(_createRayMap_sink)
        createSliceMap_sink   = Hamu.Algorithm(_createSliceMap_sink)
    else:
        createColDensMap   = Hamu.Algorithm(_createColDensMap)
        createRayMap       = Hamu.Algorithm(_createRayMap)
        createSliceMap     = Hamu.Algorithm(_createSliceMap)
        createSliceMapStar = Hamu.Algorithm(_createSliceMapStar)

    # Run for all sims
    dolengthscale = True
    isim = -1
    for ax in axes.flatten():
        ax.set_axis_off()
        ax.set_aspect("equal","datalim")
        ax.axis("off")
    for simname in simnames:
        isim += 1
        cmap  = linestyles.ColourMap(simname, hydro)
        sim   = Hamu.Simulation(simname)
        snap  = sim.FindAtTime(time)
        ax    = axes.flatten()[isim]
      
        if Slice:
            if wsink:
                data  = createSliceMap_sink(snap,hydro,los)      
            else:
                if starC:
                    data  = createSliceMapStar(snap,hydro,los)     
                else:
                    data  = createSliceMap(snap,hydro,los)     
        else: 
            if hydro == 'NH':
                if wsink: 
                    data  = createColDensMap_sink(snap,los,hydro,wsink)      
                else: 
                    data  = createColDensMap(snap,los,hydro,wsink)       
            else:
                if wsink: 
                    data = createRayMap_sink(snap,los,hydro,wsink)    
                else: 
                    data = createRayMap(snap,los,hydro,wsink)    
           
        im    = MakeImageHamu(data,hydro,wsink,ax,dolengthscale,cmap,
                       label = linestyles.Label(simname))
        dolengthscale = False

    # Add colour bar for hydro variable
    print "Making colour bar..."
    # ---- This cax works well with nrow = 1, may need adjusting, could be improved 
    cax  = fig.add_axes([0.2, 0.2, 0.4, 0.03])
    cbar = fig.colorbar(im,cax=cax,orientation="horizontal")
    label = hydrofuncs.hydro_label((hydro))
    if not "xH" in hydro:
        label = "log("+label+")"
    cbar.set_label(label,fontsize="medium",color="k")
    cbar.ax.tick_params(labelsize="medium",labelcolor="k")
    cbar.solids.set_edgecolor("face")
    # Save figure
    suffix = ""
    if len(name) > 0:
        suffix = "_"+name
    #wsname = Hamu.SimData.Settings.Settings()["CurrentWorkspace"]
    if Slice:
        folder = "../plots/vis/slice/"
        MakeDirs(folder)
        figname = folder+"slice_"+hydro+suffix+".pdf"
    else:
        folder = "../plots/vis/multiray/"
        MakeDirs(folder)
        figname = folder+"multiray_"+hydro+suffix+".pdf"
    print "Saving figure "+figname+"..."
    fig.subplots_adjust(hspace=0.0, wspace=0.0, 
                        left=0.20,right=1.0,
                        bottom=0.00,top=1.0)
    fig.set_size_inches(finches*ncols/0.8,finches*(nrows+1))
    fig.savefig(figname,
                pad_inches=0,
                dpi=dpi)
    # Crop out borders
    os.system("pdfcrop "+figname+" "+figname)
    print "Done!"

if __name__=="__main__":
    
    #time = 1.*tffcloud_code
    #MakeFigure(["IMF1_04","IMF1_05","IMF1_06","IMF1_08"],time,name="imf1_1tff",los='z',hydro='rho',Slice=True,wsink=False,starC=False)
    #MakeFigure(["IMF2_04","IMF2_05","IMF2_06","IMF2_08"],time,name="imf2_1tff",los='z',hydro='rho',Slice=True,wsink=False,starC=False)
    time = 0.5*tffcloud_code
    MakeFigure(["MASS_IMF2_02_V1","MASS_IMF2_07_V1"],time,name="v1_05tff",los='z',hydro='T',Slice=True,wsink=True,starC=False)
    #MakeFigure(["MASS_IMF2_05","MASS_IMF2_06","MASS_IMF2_08"],time,name="mass_imf2_04tff",los='z',hydro='IRflux',Slice=True,wsink=True,starC=False)
    #MakeFigure(["MASS_IMF2_04","MASS_IMF2_05","MASS_IMF2_06","MASS_IMF2_08"],time,name="mass_imf2_04tff",los='z',hydro='T',Slice=True,wsink=True,starC=False)
    #time = 0.5*tffcloud_code
    #MakeFigure(["MASS_IMF2_04","MASS_IMF2_05","MASS_IMF2_06","MASS_IMF2_08"],time,name="mass_imf2_05tff",los='z',hydro='rho',Slice=True,wsink=True,starC=False)
    #time = 0.2*tffcloud_code
    #MakeFigure(["MASS_IMF2_04","MASS_IMF2_05","MASS_IMF2_06","MASS_IMF2_08"],time,name="mass_imf2_02tff",los='z',hydro='rho',Slice=True,wsink=True,starC=False)
