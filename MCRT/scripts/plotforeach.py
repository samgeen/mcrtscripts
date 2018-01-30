'''
Plot individual figures for each output
Sam Geen, July 2014
'''
                                                                
import matplotlib as mpl
mpl.use('Agg')

import os, sys, Hamu
import numpy as np
import matplotlib.pyplot as plt

# pymses gumf
from pymses.analysis import sample_points

import spheremap

def plotspherefunc(snap,funcname,outname):
    '''
    Plot a bunch of sphere maps
    '''
    # Plotting functions for each linear slice
    def maxrho(samples,start,end):
        return samples["rho"][start:end].max()
    def rshell(samples,start,end):
        rhomax =  samples["rho"][start:end].max()
        imax = np.where(samples["rho"][start:end] == rhomax)
        radii = np.sum(samples.points[start:end,:]**2.0,axis=1)
        return np.sqrt(radii[imax])
    def cdens(samples,start,end):
        rhos =  samples["rho"][start:end]
        radii = np.sum(samples.points[start:end,:]**2.0,axis=1)
        lr = len(rhos)
        rhoc = 0.5*(rhos[0:lr-1]+rhos[1:lr])
        dr = np.abs(np.diff(radii))
        return np.sum((dr * rhoc)[rhoc > 1.5e3])
    funcs = {"maxrho": maxrho,
             "rshell": rshell,
             "cdens": cdens}
    labels = {"maxrho": "Max Density",
             "rshell": "Shell Radius",
             "cdens": "Column Density"}
    func = funcs[funcname]
    label = labels[funcname]
    amr = snap.amr_source(["rho"])
    plotlog = True
    if plotlog:
        label = "log("+label+")"
    image, imax, imin = spheremap.SphereMap(amr,func,plotlog=plotlog,radius=0.5)
    spheremap.MakeCBar(image,(imin,imax),label)
    plt.savefig("../plots/spheres/"+outname+".png",bbox_inches="tight")
    
def plotspheres(snap,simname):
    funcs = ["maxrho","rshell","cdens"]
    plotspherefuncHamu = Hamu.Algorithm(plotspherefunc)
    for func in funcs:
        plt.clf()
        os.system("mkdir -p ../plots/spheres/"+simname)
        outname = simname+"/"+func+str(snap.OutputNumber()).zfill(5)
        plotspherefuncHamu(snap,func,outname)

def run():
    Myr = 3.15569e13
    plt.clf()
    ns = ["49","48","47","00"]
    bs = ["02","00"]
    for b in bs:
        for n in ns:
            simname = "MCRTN"+n+"B"+b
            if "49" in simname:
                simname += "T3"
            if "48" in simname:
                simname += "T6"
            sim = Hamu.Simulation(simname)
            for snap in sim.Snapshots():
                plotspheres(snap,sim.Name())

if __name__=="__main__":
    run()
