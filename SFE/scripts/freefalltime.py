'''
Get the free-fall time for a simulation
Sam Geen, June 2016
'''

from startup import *

def FindVar(lines,varname):
    for line in lines:
        if varname+"=" in line or varname+" " in line:
            num = line[line.find("=")+1:]
            if num[0] == " ":
                num = num[1:]
            try:
                num = num[:num.find(" ")]
            except:
                pass
            num = float(num.replace("d","e"))
            return num    
    print "No variable",varname,"found in FindTffInNamelist"
    raise ValueError

def FindTffInNamelist(namelist,Myr):
    lines = open(namelist,"r").read().split("\n")
    tff = FindVar(lines,"trelax")
    unit = 1.0
    if Myr:
        unit = FindVar(lines,"units_time") / Myrins
    return tff * unit


def Tff(sim,Myr=False):
    snap = sim.Snapshots()[0]
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
        tff = FindTffInNamelist(namelist,Myr)
        return tff
    else:
        print "No namelist found!"
        raise OSError
