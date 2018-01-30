'''
Set the labels of each simulation in Hamu
Sam Geen, June 2015
'''

import Hamu

def run(simname):
    label = ""
    # Flux
    flux = lambda f: "S$_{*} = "+f+"$ s$^{-1}$"
    if "N00" in simname:
        label += "No photons" #flux("0") 
    else:
        label += flux("10^{"+simname[1:3]+"}")
    label += " | "
    # Supernova?
    if "NSN" in simname:
        label += "No SN"
    elif "MSN" in simname:
        label += "SN every 0.1 Myr"
    elif "HN" in simname:
        label += "Hypernova"
    elif "SN" in simname:
        label += "SN"
    # Compact?
    if "_C2" in simname:
        label += " (more compact)"
    elif "_C" in simname:
        label += " (most compact)"
    # B-field?
    if "B00" in simname:
        label += " (B=0)"
    # Delayed?
    if "_F2" in simname:
        label += " (starts at 2 t$_{ff}$)"
    if "_F3" in simname:
        label += " (starts at 3 t$_{ff}$)"
    print simname, "-->", label
    try:
        sim = Hamu.Simulation(simname)
        sim.Label(label)
    except:
        "Error setting label for simulation", simname

def runforsims(simnames):
    for simname in simnames:
        run(simname)

if __name__=="__main__":
    simnames = ["N00-NSN","N00-SN",
                "N49-NSN","N49-SN",
                "N50-NSN","N50-SN",
                "N51-NSN","N51-SN",
                "N50-MSN","N50-HN"]
    runforsims(simnames)


