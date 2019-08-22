'''
Read cubes from simulation outputs
Sam Geen, September 2018
'''

import matplotlib.pyplot as plt
import numpy as np

def readfile(hydro,nout,unit):
    '''
    Read a single cube file
    '''
    filename = hydro+"/cube_"+hydro+"_"+str(nout).zfill(5)+".npy"
    cube = np.load(filename)
    print "Read cube", hydro, "for simulation snapshot", nout
    print "Min/max =", cube.min(), cube.max(), unit
    return cube

def readcubes(nout):
    '''
    Read all the cubes and put them into dictionaries
    '''
    # List of hydrodynamic variables used
    # Bmag = sqrt(B^2), B[xyz] = B along each coordinate axis (microGauss)
    # rho = density in cm^3
    # T = temperature in K
    # xHII = hydrogen ionisation fraction from 0 to 1
    # v[xyz] = gas velocity in km/s
    # [xyz] = cartesian coordinate of cell in pc from centre of cube
    hydros = ["Bmag","Bx","By","Bz",
              "rho","T","xHII",
              "vx","vy","vz",
              "x","y","z"]
    # List of units for reference, labelling etc
    uG = "$\mu$G"
    units = [uG,uG,uG,uG,
             "cm$^{-3}$","K","(fraction)",
             "km/s","km/s","km/s",
             "pc","pc","pc"]
    # Put units and cubes into dictionaries
    # cubes["rho"] gives the cube of density, for example
    unitsdict = dict(zip(hydros, units))
    cubes = {}
    for hydro in hydros:
        cubes[hydro] = readfile(hydro,nout,unitsdict[hydro])
    return cubes, unitsdict

def plotcubes(nout):
    '''
    Make diagnostic images of the cubes
    '''
    # Read cubes
    cubes, units = readcubes(nout)
    hydros = cubes.keys()
    # Run through each hydro variable
    for hydro in hydros:
        print "Plotting image for ", hydro
        # Get the cube
        cube = cubes[hydro]
        # Flatten cube to find maximum along line of sight
        maxcube = np.max(cube,axis=1)
        if maxcube.min() > 0.0:
            maxcube = np.log10(maxcube)
            label = "log(max "+hydro+")"
        else:
            label = "max "+hydro
        # Plot
        plt.clf()
        plt.imshow(maxcube,cmap = "RdPu_r")
        plt.colorbar(label=label+" / "+units[hydro])
        plt.savefig("images/maximage_"+hydro+"_"+str(nout).zfill(5)+".pdf")
        
if __name__=="__main__":
    # Only made cubes for snapshot 38 so far
    #cubes, units = readcubes(50)
    plotcubes(39)
    
