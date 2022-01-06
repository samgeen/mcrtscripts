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

cmap = option_d.test_cm
import sys
sys.path.append("../../scripts/")
import hydrofuncs
#cmap = brewer2mpl.get_map('RdPu', 'Sequential', 9,reverse=True).mpl_colormap
#cmap = brewer2mpl.get_map('OrRd', 'Sequential', 9,reverse=True).mpl_colormap

<<<<<<< HEAD
=======
isims = {}
isim = 0

def reset():
    global isims, isim
    isims = {}
    isim = 0

>>>>>>> dacd3981783136068743eff48673f6efc6d88fab
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
<<<<<<< HEAD
=======
    global isims, isim
>>>>>>> dacd3981783136068743eff48673f6efc6d88fab
    colours = ["#a6cee3",
 	       	"#1f78b4",
		"#b2df8a",
		"#33a02c",
		"#fb9a99",
		"#e31a1c",
		"#fdbf6f",
		"#ff7f00",
                "#ff5900"
		]
<<<<<<< HEAD
    try:
        return colours[allsims.index(simname)]
=======
    if not simname in isims:
        isims[simname] = isim
        isim += 1
    try:
        return colours[isims[simname]]
>>>>>>> dacd3981783136068743eff48673f6efc6d88fab
    except:
        print("No colour for simulation", simname,"found, returning #000000")
        return "#000000"

def Linestyle(simname):
    lines = ["--","-",":","-.","-","-","--"]
    return "-"

def Label(simname):
    label = simname.replace("_35MSUN","")
    label = label.replace("_","\_") # TODO: fix this
    return label
