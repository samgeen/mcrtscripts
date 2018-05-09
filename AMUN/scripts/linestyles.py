'''
Set up the line styles for the simulations
Sam Geen, December 2017
'''

import customplot
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from collections import OrderedDict


import prettyplotlib as ppl
from prettyplotlib import plt
from prettyplotlib import brewer2mpl

#cmap = brewer2mpl.get_map('RdPu', 'Sequential', 9,reverse=True).mpl_colormap
cmap = brewer2mpl.get_map('OrRd', 'Sequential', 9,reverse=True).mpl_colormap

# SIMULATIONS
# 01 - NOUV/NOWINDS
# 02 -   UV/NOWINDS
# 03 - NOUV/  WINDS
# 04 -   UV/  WINDS

def ColourMap(simname=None):
    return cmap

def RunNum(simname):
    run, num = simname.split("_")
    num = int(num)
    return run, num

def Colour(simname):
    # Colour palette is 4-class Paired from ColorBrewer2.org
    colours = ["#a6cee3",
               "#1f78b4",
               "#b2df8a",
               "#33a02c"]
    run, num = RunNum(simname)
    return colours[num-1]

def Linestyle(simname):
    run, num = RunNum(simname)

def Label(simname):
    run, num = RunNum(simname)
    labels = ["No Feedback",
              "UV Only",
              "Winds Only",
              "Winds \& UV"]
    return labels[num-1]
