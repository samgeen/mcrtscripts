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

import option_d
import option_c

CURRSIMSET = None

cmap = option_d.test_cm
import sys
sys.path.append("../../scripts/")
import hydrofuncs
#cmap = brewer2mpl.get_map('RdPu', 'Sequential', 9,reverse=True).mpl_colormap
#cmap = brewer2mpl.get_map('OrRd', 'Sequential', 9,reverse=True).mpl_colormap

simlabels = {}

isims = {}
isim = 0

def reset():
    global isims, isim
    isims = {}
    isim = 0

def ColourMap(simname=None, hydro='rho'):
    return hydrofuncs.cmap(hydro)
    '''
    if hydro == 'rho':
        cmap = option_d.test_cm
    elif hydro == 'nH':
        cmap = option_d.test_cm
    elif hydro == 'NH':
        cmap = option_d.test_cm
    elif hydro == 'T':
        #cmap = option_c.test_cm
        cmap = "jet" 
    else:
        cmap = brewer2mpl.get_map('OrRd', 'Sequential', 9,reverse=True).mpl_colormap
    return cmap
    '''

def Colour(simname):
    # Colour palette is 4-class Paired from ColorBrewer2.org
    #colours = ["#a6cee3",
    #           "#1f78b4",
    #           "#b2df8a",
    #           "#33a02c"
    #          ]
    global isims, isim
    colours = ["#a6cee3",
               "#1f78b4",
               "#b2df8a",
               "#33a02c",
               "#fb9a99",
               "#e31a1c",
               "#fdbf6f",
               "#ff7f00",
               "#ff5900"]
    if not simname in isims:
        isims[simname] = isim
        isim += 1
    try:
        return colours[isims[simname]]
    except:
        print("No colour for simulation", simname,"found, returning #000000")
        return "#000000"

def Linestyle(simname):
    lines = ["--","-",":","-.","-","-","--"]
    return "-"

def Label(simname):
    global CURRSIMSET
    if CURRSIMSET is not None:
        label = simlabels[CURRSIMSET][simname]
    else:
        label = simname.replace("_35MSUN","")
        label = label.replace("_","\_") # TODO: fix this
    return label
