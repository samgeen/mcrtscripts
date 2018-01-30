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

icmap = brewer2mpl.get_map('RdPu', 'Sequential', 9,reverse=True).mpl_colormap
imfmap = brewer2mpl.get_map('OrRd', 'Sequential', 9,reverse=True).mpl_colormap

def Colour(simname):
    if "IMF" in simname:
        colmap = imfmap
    else:
        colmap = icmap
    num = float(simname[-2:])
    # TODO: REMOVE TOO-LIGHT COLOURS
    col = colmap(num/13.0)
    return col

