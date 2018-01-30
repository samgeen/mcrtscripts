'''
Uniform line styles and colours for the project
Sam Geen, January 2016
'''

import Hamu
import customplot
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from collections import OrderedDict

'''
Sims:
Winds? - NW, W
RT? - RT, NRT
'''

def finditem(tosearch,simname):
    for key in tosearch.iterkeys():
        if key == simname:
            return tosearch[key]
    # Oops! Not found the result?
    print "Error: No result for", simname, "in keys",tosearch.keys()
    raise KeyError

colours = {}
winds = ["W","NW"]
# Pink colour scheme
cols = filter(None,'''
#fbb4b9
#f768a1
#c51b8a
#7a0177
'''.split("\n"))
# Varied colours with different lightness scheme
cols = ['#000000', '#2a76d4']
colours = OrderedDict()
for key, value in zip(winds,cols):
    colours[key] = value

def colour(simname):
    return finditem(colours,simname.split("-")[0])

lines = {}
lines["RT"] = "-"
lines["NRT"] = "--"

def line(simname):
    return finditem(lines,simname.split("-")[1])

windlabels = {}
#sizelabels["L"] = r"$R/R_0 = 2.25$"
#sizelabels["M"] = r"$R/R_0 = 1$"
#sizelabels["S"] = r"$R/R_0 = 0.5625$"
#sizelabels["XS"] = r"$R/R_0 = 0.25$"
windlabels["NW"] = r"No winds"
windlabels["W"]  = r"Winds"

def windlabel(simname):
    return finditem(windlabels,simname.split("-")[0])

rtlabels = {}
rtlabels["RT"]  = r"Radiation Included"
rtlabels["NRT"] = r"No radiation"

def rtlabel(simname):
    return finditem(rtlabels,simname.split("-")[0])

def windlegend():
    plotlines = []
    for size, col in colours.iteritems():
        plotlines.append(mlines.Line2D([], [], 
                                       color=col,
                                       linestyle="-",
                                       label=windlabel(size)))
    labels = [line.get_label() for line in plotlines]
    return plotlines, labels

def rtlegend():
    plotlines = []
    for rt, line in lines.iteritems():
        plotlines.append(mlines.Line2D([], [], 
                                       color="k",
                                       linestyle=line,
                                       label=rtlabel(rt)))
    labels = [line.get_label() for line in plotlines]
    return plotlines, labels


