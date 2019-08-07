'''
Make a data cube from the fields in a simulation
Sam Geen, May 2015
'''

import copy

from startup import *

from astropy.io import fits

import stellars

def run(snap,folder,hydros,pos,radius):
    # Open snapshot
    ro = snap.RawData()
    print "Sampling grid ...",
    hydros = copy.copy(hydros)
    amr = hydrofuncs.amr_source(ro,hydros)
    # Make grid
    lmin = ro.info["levelmin"]
    lsize = 128
    boxlen = ro.info["boxlen"]
    pos = np.array(pos)/boxlen
    coords = np.linspace(-0.5,0.5,lsize)*2.0*radius/boxlen
    grid = np.meshgrid(coords+pos[0],coords+pos[1],coords+pos[2])
    points = np.array([grid[0].flatten(),grid[1].flatten(),grid[2].flatten()])
    points = points.T
    print points.min(), points.max(), points.shape
    samples = pymses.analysis.sample_points(amr,points)
    # Sample for each hydro variable
    hydros += ["x","y","z"]
    for hydro in hydros:
        filename = folder+"/"+hydro+\
                   "/cube_"+hydro+"_"+str(snap.OutputNumber()).zfill(5)+".fits"
        if not os.path.exists(filename):
            print "Making", filename,"...",
            if hydro not in "xyz":
                scale = hydrofuncs.scale_by_units(ro,hydro)
                hydrocube = scale(samples)
                hydrocube = hydrocube.reshape([lsize,lsize,lsize])
            else:
                if hydro == "x":
                    hydrocube = (grid[0]-pos[0])
                if hydro == "y":
                    hydrocube = (grid[1]-pos[1])
                if hydro == "z":
                    hydrocube = (grid[2]-pos[2])
                hydrocube *= boxlen
            print "Minmax=",hydrocube.min(),hydrocube.max(),"..."
            #np.save(filename, hydrocube.astype("float32"))
            hdu = fits.PrimaryHDU(hydrocube)
            hdul = fits.HDUList([hdu])
            hdul.writeto(filename)
        else:
            print filename, "exists, ignoring"
    print "Done!"

def makedir(folder):
    try:
        os.mkdir(folder)
    except:
        pass

def runforsim(sim,nouts=None,times=None,pos=None,radius=None):
    hydros = ["vx","vy","vz","rho","T","xHII","NpIR","NpOpt","NpHII","NpHeII","NpHeIII","P","Pnontherm","xHeII"]
    simname = sim.Name()
    simfolder = "../cubes/fits2/"+simname
    makedir(simfolder)
    snap = sim.Snapshots()[0]
    myr = snap.RawData().info["unit_time"].express(C.Myr)
    if nouts is not None:
        times = []
        timedict = {snap.OutputNumber():snap.Time() for snap in sim.Snapshots()}
        for nout in nouts:
            times.append(timedict[nout])
    else:
        times = [t/myr for t in times]
    for hydro in hydros:
        makedir(simfolder+"/"+hydro)
    for cpos in "xyz":
        makedir(simfolder+"/"+cpos)
    if times is None:
        for snap in sim.Snapshots():  
            run(snap,folder=simfolder,hydros=hydros,pos=pos,radius=radius)
    else:
        for time in times:
            snap = sim.FindAtTime(time)
            run(snap,folder=simfolder,hydros=hydros,pos=pos,radius=radius)


if __name__=="__main__":
    # Use IMF2, winds + UV
    sim = hamusims["ERIC-WIND"]
    # Pick the last output - TODO, CHANGE TO SPECIFIC OUTPUT!!!
    #nout = snap.OutputNumber()
    #nouts = [nout]
    #snap1 = sim.Snapshots()[17]
    snap2 = sim.Snapshots()[-1]
    myr = snap2.RawData().info["unit_time"].express(C.Myr)
    times = [snap.Time()*myr for snap in [snap2]]
    #times = np.linspace(0.0,1.0,11)+3.38 # Myr
    #snap = sim.FindAtTime(times[0]/myr)
    # Pick the first star
    stars = stellars.FindStellar(snap)
    pos = [stars.x[0],stars.y[0],stars.z[0]]
    radius = 25.0
    #pos = np.zeros(3)+0.5
    #radius = 0.4
    runforsim(sim,times=times,pos=pos,radius=radius)
