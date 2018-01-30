'''
Make a data cube from the fields in a simulation
Sam Geen, May 2015
'''

import Hamu
import numpy as np
import pymses
import hydrofuncs

import os

from pymses.utils import constants as C

def run(snap,folder,hydro):
    ro = snap.RawData()
    filename = folder+"/"+hydro+\
        "/cube_"+hydro+"_"+str(snap.OutputNumber()).zfill(5)+".npy"
    if os.path.exists(filename):
        print filename, "exists, returning"
        return
    else:
        print "Making", filename,"...", 
    lmin = ro.info["levelmin"]
    lsize = 2**lmin
    coords = np.arange(0.0,1.0,1.0/float(lsize))
    grid = np.meshgrid(coords,coords,coords)
    points = np.array([grid[0].flatten(),grid[1].flatten(),grid[2].flatten()])
    points = points.T
    print "Sampling grid ...", 
    amr = hydrofuncs.amr_source(ro,hydro)
    scale = hydrofuncs.scale_by_units(ro,hydro)
    samples = pymses.analysis.sample_points(amr,points)
    dens = scale(samples)
    np.save(filename, dens.astype("float32"))
    print "Done!"

def makedir(folder):
    try:
        os.mkdir(folder)
    except:
        pass

def runforsim(simname,times=None):
    hydros = ["vx","vy","vz","rho","T","xHII","Bmag"]
    sim = Hamu.Simulation(simname)
    simfolder = "../cubes/"+sim.Name()
    makedir(simfolder)
    myr = None
    for hydro in hydros:
        makedir(simfolder+"/"+hydro)
    if times is None:
        for snap in sim.Snapshots():  
            for hydro in hydros:
                run(snap,folder=simfolder,hydro=hydro)
    else:
        for time in times:
            if myr is None:
                snap = sim.Snapshots()[0]
                myr = snap.RawData().info["unit_time"].express(C.Myr)
            snap = sim.FindAtTime(time/myr)
            for hydro in hydros:
                run(snap,folder=simfolder,hydro=hydro)


if __name__=="__main__":
    #runforsim("N48_M4_B02_C2")
    Hamu.Workspace("HIISN")
    tff = 2.53 + 3.0
    times = tff + np.arange(0,1,0.25)
    runforsim("N50-SN", times)
