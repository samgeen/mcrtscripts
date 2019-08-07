'''
Make a data cube from the fields in a simulation in FITS format with yt
Sam Geen, July 2019
'''

from startup import *

from astropy.io import fits
import yt

import stellars

def run(snap,folder):
    ro = None
    for hydro in ["density","x-acceleration","y-acceleration","z-acceleration"]:
        filename = folder+"/"+hydro+\
                   "/cube_"+hydro+"_"+str(snap.OutputNumber()).zfill(5)+".fits" 
        if not os.path.exists(filename):
            if ro is None:
                # Open snapshot
                ro = snap.RawData()
                boxlen = ro.info["boxlen"]
                utime = snap.RawData().info["unit_time"].express(C.s)
                ulength = snap.RawData().info["unit_length"].express(C.cm)
                uacc = ulength / utime**2.0
                # This is kind of hacky, but use Pymses snapshot properties to load yt
                numstr = str(ro.iout).zfill(5)
                path = ro.output_repos+"output_"+numstr+"/info_"+numstr+".txt"
                #fields = ["Density",
                #          ("gravity","x-acceleration"),
                #          ("gravity","y-acceleration"),
                #          ("gravity","z-acceleration")]
                ds = yt.load(path)#,fields=fields)
                cube = ds.covering_grid(level=2, left_edge=[0,0.0,0.0],
                                        dims=ds.domain_dimensions*2**2)

            if "acceleration" in hydro:
                hydrocube = cube[hydro]*uacc
                hydrocube /= boxlen # Hack fix to issue with yt
            else:
                hydrocube = cube[hydro]
            makedir(folder+"/"+hydro)
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

def runforsim(sim,nouts=None,times=None):
    simname = sim.Name()
    simfolder = "../cubes/fits/"+simname
    makedir(simfolder)
    snap = sim.Snapshots()[0]
    myr = snap.RawData().info["unit_time"].express(C.Myr)
    if nouts is not None:
        times = []
        timedict = {snap.OutputNumber():snap.Time() for snap in sim.Snapshots()}
        for nout in nouts:
            times.append(timedict[nout])
    elif times is not None:
        times = [t/myr for t in times]
    if times is None:
        for snap in sim.Snapshots():  
            run(snap,folder=simfolder)
    else:
        for time in times:
            snap = sim.FindAtTime(time)
            run(snap,folder=simfolder)


if __name__=="__main__":
    # Use IMF2, winds + UV
    sim = hamusims["IC06"]
    # Pick the last output - TODO, CHANGE TO SPECIFIC OUTPUT!!!
    #nout = snap.OutputNumber()
    #nouts = [nout]
    #snap = sim.Snapshots()[-1]
    #myr = snap.RawData().info["unit_time"].express(C.Myr)
    #times = np.linspace(0.0,1.0,11)+3.38 # Myr
    #snap = sim.FindAtTime(times[0]/myr)
    # Pick the first star
    #stars = stellars.FindStellar(snap)
    #pos = [stars.x[0],stars.y[0],stars.z[0]]
    #radius = 25.0
    for sim in hamusims.itervalues():
        runforsim(sim)
