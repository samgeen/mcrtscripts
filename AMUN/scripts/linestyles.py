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

def _RunIndex(simname):
    d = {""}

def Colour(simname):
    # Colour palette is 4-class Paired from ColorBrewer2.org
    colors = ["#a6cee3",
              "#1f78b4",
              "#b2df8a",
              "#33a02c"]
    run, num = RunNum(simname)
    if run == "imf":
        colmap = imfmap
    else:
        colmap = icmap
    # TODO: REMOVE TOO-LIGHT COLOURS
    col = colmap(float(num)/13.0)
    return col

def Label(simname,sub=False):
    run, num = RunNum(simname)
    if run == "imf":
        label = imfabbrevs[num-1]
        if sub:
            label += starsub
        return label
    else:
        label = icabbrevs[num-1]
        if sub:
            label += turbsub
        return label
    print "Oops, something went wrong in linestyles.Labels for simname",simname
    raise ValueError

#def Legend(run,scatter=True):
#    if run == "imf":
        
