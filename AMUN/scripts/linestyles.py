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
#cmap = brewer2mpl.get_map('RdPu', 'Sequential', 9,reverse=True).mpl_colormap
#cmap = brewer2mpl.get_map('OrRd', 'Sequential', 9,reverse=True).mpl_colormap

# SIMULATIONS
# 01 - NOUV/NOWINDS
# 02 -   UV/NOWINDS
# 03 - NOUV/  WINDS
# 04 -   UV/  WINDS
# 05 -   UV+press/  WINDS
# 06 -   UV+IR+press/  WINDS

def ColourMap(simname=None, hydro='rho'):
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

def RunNum(simname):
    if ('V' in simname):
        if ('02' in simname):
            run = '02'
            num = 2 
        elif ('07' in simname):
            run = '07'
            num = 7 
        elif ('09' in simname):
            run = '09'
            num = 9 
    else:
        run, num = simname.split("_")[-2], simname.split("_")[-1]
    num = int(num)
    return run, num

def Colour(simname):
    # Colour palette is 4-class Paired from ColorBrewer2.org
    #colours = ["#a6cee3",
    #           "#1f78b4",
    #           "#b2df8a",
    #           "#33a02c"
    #          ]
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
    run, num = RunNum(simname)
    return colours[num-1]

def Linestyle(simname):
    run, num = RunNum(simname)

def Label(simname):
    run, num = RunNum(simname)
    labels = ["No Feedback",  # 01
              "UV",           # 02
              "Winds",        # 03
              "Winds + UV",   # 04
              "Winds + UVp",  # 05
              "Winds + UVp + IRp", # 06
              "UVp",          # 07
              "UVp + IRp",    # 08
              "UVp + Boost"]  # 09
    return labels[num-1]
