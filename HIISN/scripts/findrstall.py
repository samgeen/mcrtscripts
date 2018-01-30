'''
Find rstall in the HIISN model
Sam Geen, July 2015
'''

import sys
sys.path.append("/home/sgeen/MC_RT/scripts")

import numpy as np

import Hamu
Hamu.Workspace("HIISN")

import outflowmodel, ragacompare, profilesphere, testprofilefit
import edgeradius

flux = 0.0

class Prof(object):
    def __init__(self,snap):
        #self._n0 = outflowmodel.FindnH(sim)
        #rhofit,w = testprofilefit.Run(sim.Name(),outflowmodel.tstart)
        #r0 = np.exp((np.log(rhofit) - np.log(self._n0))/w)
        # HACK - HAND-FITTED VALUES
        w = 0.0 # INSIDE CLOUD CORE
        r0 = 3.6 # HACK 5.0 # based on testprofiel
        self._n0 = 1612.0 #HACK 612.0
        r0 *= outflowmodel.pcincm
        self._r0 = r0
        self._r0w = r0**w
        self._minusw = -w

def edgerange(snap):
    edges, denses = edgeradius.FindEdgesHamu(snap)
    elims = [edges[edges > 0.0].min(), edges.max()]
    return elims

def stallrad(snap):
    prof = Prof(snap)
    w = -prof._minusw # eh
    r0w = prof._r0w
    r0 = prof._r0
    ci = ragacompare.cs
    return anarstall(flux,prof._n0,r0,w,ci)

def anarstall(S,n0,r0,w,ci):
    G = ragacompare.G
    mH = ragacompare.mH
    X = ragacompare.X
    rho0 = n0 * mH / X
    r0w = r0**w
    alphaB = ragacompare.beta2
    flux = S
    #vfact = 2.0 # Escape velocity
    vfact = 6.0/5.0 # Virial velocity
    rstall = ((3.0*flux / (4*np.pi * alphaB))**0.5 * \
                  (1.0/r0w**2) * (3 - w) * (ci**2) *mH/X / \
                  (4*vfact * np.pi * G * rho0**2))**(2.0 / (7.0 - 4*w))
    #import pdb; pdb.set_trace()
    return rstall/outflowmodel.pcincm

def run(simname):
    sim = Hamu.Simulation(simname)
    snap = sim.Snapshots()[5] # Start of source
    outflowmodel.tstart = snap.Time()
    rstall = stallrad(snap)
    prof = Prof(snap)
    r0 = prof._r0/outflowmodel.pcincm
    boxlen = snap.RawData().info["boxlen"]
    print "Stall radius @ " ,flux, " photons/s: ", rstall, "pc"
    print "Box length: ", boxlen
    er = edgerange(snap)
    print "Rough cloud radius:", er
    print "rstall / rcloud_min:", rstall / er[0]
    #print "rstall / rcloud_min * (0.4/0.5)**2:", rstall / er[0] * (4.0/5.0)**2
    print "rstall / r0:",rstall / r0
    print "----------"

if __name__=="__main__":
    simname = "N00-NSN"
    flux = 1e49
    run(simname)
    flux = 1e50
    run(simname)
    flux = 1e51
    run(simname)
