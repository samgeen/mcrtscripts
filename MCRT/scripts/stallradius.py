'''
Calculate the stall radius (as well as by a function of cloud radius)
Sam Geen, June 2015
'''

import numpy as np
import ragacompare, outflowmodel, Hamu

def stallradius(sim):
    prof = ragacompare.ProfAna(sim)
    w = -prof._minusw # eh
    r0w = prof._r0w
    mH = ragacompare.mH
    X = ragacompare.X
    rho0 = prof._n0 * mH / X
    alphaB = ragacompare.beta2
    G = ragacompare.G
    ci = ragacompare.cs
    flux = outflowmodel.FindFlux(sim)
    rstall = ((3.0*flux / (4*np.pi * alphaB))**0.5 * \
                  (1.0/prof._r0w**2) * (3 - w) * ci**2 *mH/X / \
                  (8 * np.pi * G * rho0**2))**(2.0 / (7.0 - 4*w))
    #import pdb; pdb.set_trace()
    return rstall/outflowmodel.pcincm

def rcloud(simname,limit="average"):
    rmin = 3.0
    rmax = 12.0
    if limit == "average":
        rcloud = 0.5*(rmin+rmax)
    elif limit == "min":
        rcloud = rmin
    elif limit == "max":
        rcloud = rmax
    else:
        print "Wrong limit!"
        raise ValueError
    if "_C2" in simname:
        rcloud *= 0.75**2
    elif "_C" in simname:
        rcloud *= 0.5**2
    return rcloud

def run(simnames):
    sims = {}
    rstalls = {}
    rclouds = {}
    for simname in simnames:
        sims[simname] = Hamu.Simulation(simname)
        rclouds[simname] = rcloud(simname,limit="min")
        rstalls[simname] = stallradius(sims[simname])
    print "RCLOUD, RSTALL, RSTALL / RCLOUD"
    print "--------------------"
    for simname in simnames:
        print simname, rclouds[simname], rstalls[simname],\
            rstalls[simname] / rclouds[simname]

if __name__=="__main__":
    simnames = ["N48_M4_B02","N47_M4_B02","N49_M4_B02",
                "N48_M4_B02_C2","N48_M4_B02_C"]
    run(simnames)
