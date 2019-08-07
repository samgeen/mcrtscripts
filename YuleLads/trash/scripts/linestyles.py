'''
Uniform line styles and colours for the project
Sam Geen, January 2016
'''

#from startup import *

import numpy as np
#import customplot
import simlabels

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from collections import OrderedDict

'''
Sims:
1-13
'''

allsims = simlabels.allsims

def finditem(tosearch,simname):
    for key in tosearch.iterkeys():
        if key == simname:
            return tosearch[key]
    # Oops! Not found the result?
    print "Error: No result for", simname, "in keys",tosearch.keys()
    raise KeyError

colours = {}
# Varied colours with different lightness scheme

reds = np.arange(0,1.01,1.0/float(len(allsims)-1))
cols = [(r,0.0,0.0) for r in reds]
colours = OrderedDict()
for key, value in zip(allsims,cols):
    colours[key] = value

def colour(simname):
    return finditem(colours,simname)

# Linestyles not really used yet !!!

lines = {}
lines["UVP"] = "-"
lines["UVNP"] = "--"
lines["NUV"] = ":"

def line(simname):
    return "-"
#    return finditem(lines,simname.split("-")[1])

#irlabels = {}
#irlabels["IRP"] = r"IR Pressure"
#irlabels["IRNP"]  = r"IR, No Pressure"
#irlabels["NIR"]  = r"No IR"

#def irlabel(simname):
#    return finditem(irlabels,simname.split("-")[0])

# EDIT BELOW HERE

uvlabels = {}
uvlabels["UVP"]  = r"UV Photoionisation + Pressure"
uvlabels["UVNP"] = r"UV Photoionisation only"
uvlabels["NUV"] = r"No UV"

def uvlabel(simname):
    return finditem(uvlabels,simname.split("-")[0])

def simlegend():
    plotlines = []
    for sim, col in colours.iteritems():
        plotlines.append(mlines.Line2D([], [], 
                                       color=col,
                                       linestyle="-",
                                       label=sim))
    labels = [line.get_label() for line in plotlines]
    return plotlines, labels

def uvlegend():
    plotlines = []
    for rt, line in lines.iteritems():
        plotlines.append(mlines.Line2D([], [], 
                                       color="k",
                                       linestyle=line,
                                       label=uvlabel(rt)))
    labels = [line.get_label() for line in plotlines]
    return plotlines, labels


