'''
Print power law profile stuff
Sam Geen, April 2015
'''

import Hamu, outflowmodel

pcincm = 3.08567758e18

def RunForSim(sim):
    r0 = outflowmodel.FindScaleRadius(sim)
    A = outflowmodel.FindProfilePower(sim)
    rho0 = outflowmodel.FindnH(sim,rcut=0.1)
    return r0/pcincm, A, rho0

def Run(simnames):
    sims = []
    params = []
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        sims.append(sim)
        params.append(RunForSim(sim))

    print "OKAY LET'S GO NOW YES"
    print ""
    print "YEAH"
    print ""

    for param, sim in zip(params,sims):
        print sim.Name(), "r0, A, rho0", param


if __name__=="__main__":
    simnames = ["N48_M4_B02",
                "N48_M4_B00",
                "N48_M4_B02_C",
                "N48_M4_B02_C2",
                "N48_M4_B02_F2",
                "N48_M4_B02_F3"]
    Run(simnames)
