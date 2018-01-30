'''
Uniform line styles and colours for the project
Sam Geen, January 2016
'''

import Hamu
import customplot
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from collections import OrderedDict

sims = ["N00-NSN","N00-SN",
        "N49-NSN","N49-SN",
        "N50-NSN","N50-SN",
        "N51-NSN","N51-SN",
        "N50-MSN","N50-HN"]

# NSN = dashed
# SN = solid
# HN = dot
# MSN = dash-dot

# Photon emission rates = orange-red
# Colour table YlOrRd

'''
EMISSION RATE PROPERTIES
'''

# Make list of emissions
def makeems():
    ems = OrderedDict()
    ws = Hamu.CurrentWorkspace()
    if ws == "HIISN":
        # Old values 
        #ems["N51"] = "#bd0026"
        #ems["N50"] = "#f03b20"
        #ems["N49"] = "#fd8d3c"
        #ems["N00"] = "#fecc5c"
        ems["N51"] = "#feb24c"
        ems["N50"] = "#e31a1c"
        ems["N49"] = "#800026"
        ems["N00"] = "#111111"
    elif ws == "HIISN4":
        # TODO: CHANGE THIS TO ANOTHER COLOUR SCHEME!!!
        ems["N49"] = "#bd0026"
        ems["N48"] = "#f03b20"
        ems["N47"] = "#fd8d3c"
        ems["N00"] = "#fecc5c"
    return ems

def emlabel(simname):
    if "00" in simname:
        return "No photons"
    else:
        return "$10^{"+simname[1:3]+"}$ /s"

def emlegend():
    ems = makeems()
    plotlines = []
    for em, col in ems.iteritems():
        plotlines.append(mlines.Line2D([], [], 
                                       color=col,
                                       linestyle="-",
                                       label=emlabel(em)))
    labels = [line.get_label() for line in plotlines]
    return plotlines, labels

# Convenience functions
def col(simname):
    ems = makeems()
    # Match the simulation name to a colour
    for key, val in ems.iteritems():
        if key in simname:
            return val
    # A default value
    return "#000000"

'''
SN PROPERTIES
'''

# Make list of SN properties
def makesns():
    sns = OrderedDict()
    ws = Hamu.CurrentWorkspace()
    if ws == "HIISN":
        sns["-NSN"] = "--"
        sns["-SN" ] = "-"
        sns["-HN" ] = ":"
        sns["-MSN"] = "-."
        sns["-NE"] = ":" # Figure out a better way for this
        sns["-ME"] = "--" # Figure out a better way for this
        sns["-MNE"] = "-." # Figure out a better way for this
    elif ws == "HIISN4":
        sns["-NSN"] = "--"
        sns["-SN" ] = "-"
    return sns

def snlabel(simname):
    label = ""
    if "-NSN" in simname:
        label += "No SN"
    elif "-MSN" in simname:
        label += "SN every 0.1 Myr"
    elif "-HN" in simname:
        label += "Hypernova"
    elif "-NE" in simname:
        label += "SN (no ejecta)"
    elif "-ME" in simname:
        label += "SN (momentum injection)"
    elif "-MNE" in simname:
        label += "SN (momentum, no ejecta)"
    elif "-SN" in simname:
        label += "SN"
    return label

def line(simname):
    sns = makesns()
    for key, val in sns.iteritems():
        if key in simname:
            return val
    # Default option
    return "-"

def snlegend():
    sns = makesns()
    plotlines = []
    for sn, ln in sns.iteritems():
        plotlines.append(mlines.Line2D([], [], 
                                       color="k",
                                       linestyle=ln,
                                       label=snlabel(sn)))
    labels = [line.get_label() for line in plotlines]
    return plotlines, labels
    

if __name__=="__main__":
    simnames = cols.keys()
    isim = 1
    for simname in simnames:
        print simname, col(simname), line(simname)
        plt.plot([0,1],[isim,isim],
                 color=col(simname),
                 linestyle=line(simname),
                 label=simname)
        isim += 1
    plt.ylim([0,isim])
    plt.legend(fontsize="xx-small",frameon=False)
    plt.savefig("../plots/linestyletester.pdf")
