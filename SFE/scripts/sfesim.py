'''
Find the SFE of a snapshot according to the simulation
Sam Geen, June 2016
'''

from startup import *

import sinks, sinkextinction

def FindCloudMassInNamelist(namelist):
    lines = open(namelist,"r").read().split("\n")
    for line in lines:
        if "mass_c=" in line or "mass_c " in line:
            num = line[line.find("=")+1:]
            try:
                num = num[:num.find(" ")]
            except:
                pass
            num = float(num)
            return num
    print "No cloud mass found in FindCloudMassInNamelist"
    raise ValueError

def CloudMassInSim(snap):
    cloudmass = 1e4 # DEFAULT VALUE!
    ro = snap.RawData()
    namelists = glob.glob(ro.output_repos+"/*.nml")
    # Catch potential problems reading more than one namelist
    if len(namelists) > 1:
        print "You nave more than one namelist in this folder"
        print "I can't tell which one I should read"
        print "Namelists found:", namelists
        print "Either delete/rename the ones you don't need or modify sfe.py"
        raise RuntimeError
    namelist = namelists[0]
    if os.path.exists(namelist):
        #try:
        cloudmass = FindCloudMassInNamelist(namelist)
        #except:
        #    pass
    return cloudmass
    

def SFEsim(snap):
    '''
    Find the total 'SFE' of a snapshot
    This is just the sink mass / 1e4 solar masses
    '''
    # Find cloud mass
    cloudmass = CloudMassInSim(snap.hamusnap)
    sinkmass = np.sum(sinks.FindSinks(snap.hamusnap).mass)
    # Note that the cloudmass is the *initial* mass, not the current one
    return sinkmass / cloudmass

SFEsimHamu = Hamu.Algorithm(SFEsim)
        

if __name__=="__main__":
    sim = Hamu.Simulation("L-RT")
    snap = sim.Snapshots()[-1]
    print SFEsimHamu(snap)
    print sinks.FindSinks(snap).mass
