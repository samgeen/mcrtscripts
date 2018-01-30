'''
Make profiles with and without magnetic fields at tstart
Sam Geen, April 2015
'''

import customplot
import matplotlib.pyplot as plt
import numpy as np
import Hamu

import hydrofuncs

import pymses
from pymses.utils import constants as C

import profilesphere
profilemodule = profilesphere

import linestyles

def run():
   tstart = 1.25
   simnames = ["N00_M4_B02","N48_M4_B00"]
   sims = [Hamu.Simulation(simname) for simname in simnames]
   # Rho
   plt.clf()
   for sim in sims:
       firstsnap = sim.Snapshots()[0]
       myr = firstsnap.RawData().info["unit_time"].express(C.Myr)
       snap = sim.FindAtTime(tstart/myr)
       line = "-"
       label = "B-field"
       if "B00" in sim.Name():
           label = "No B-Field"
       r, prof = profilemodule.profileHamu(snap,"rho",1e6)
       plt.plot(r,prof,color=linestyles.col(sim.Name()),
                linestyle=line,label=label)
   plt.text(1,2,"t = "+str(tstart)+" Myr")
   plt.yscale("log")
   plt.xlabel("Time / Myr")
   plt.ylabel(hydrofuncs.hydro_label("rho"))
   plt.legend(fontsize="xx-small",frameon=False)
   plt.savefig("../plots/magnomag/magnomag_rho.pdf")
   # Infall
   plt.clf()
   for hydro in ["vrad","spd"]:
       for sim in sims:
           firstsnap = sim.Snapshots()[0]
           myr = firstsnap.RawData().info["unit_time"].express(C.Myr)
           snap = sim.FindAtTime(tstart/myr)
           line = "--"
           label = "B-field"
           if "B00" in sim.Name():
               label = "No B-Field"
           r, prof = profilemodule.profileHamu(snap,hydro,1e6)
           if hydro == "vrad":
               line = "-"
               plt.plot(r,prof,color=linestyles.col(sim.Name()),
                        linestyle=line,label="Radial Velocity ("+label+")")
           else:
              line = "--"
              plt.plot(r,prof,color=linestyles.col(sim.Name()),
                       linestyle=line,label="Total Velocity ("+label+")")
       #plt.plot(r,-prof,"r"+line)
   plt.yscale("linear")
   plt.xlabel("Time / Myr")
   plt.ylabel(hydrofuncs.hydro_label("vrad"))
   plt.legend(loc="lower right",fontsize="xx-small",ncol=2,frameon=False)
   plt.savefig("../plots/magnomag/magnomag_vrad.pdf")

if __name__=="__main__":
    run()
