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

imfabbrevs = "STE,GIL,STU,THV,POT,ASK,HUR,SKY,BJU,GLU,GAT,KET,KER".split(",")
icabbrevs = "GRY,JOL,JOU,GAL,TIO,CAG,SNE,BEF,STL,MAR,TAN,OLD,YAG".split(",")
starsub = r"$_{\normalfont\textsc{stars}}$"
turbsub = r"$_{\normalfont\textsc{turb}}$"
starsc  = r"\textsc{stars}"
turbsc  = r"\textsc{turb}"

def ColourMap(simname):
    if "imf" in simname.lower():
        return imfmap
    else:
        return icmap

def RunNum(simname):
    if "IMF" in simname:
        run = "imf"
    else:
        run = "ic"
    num = 7
    return run, num

def Colour(simname):
    run, num = RunNum(simname)
    if run == "imf":
        colmap = imfmap
    else:
        colmap = icmap
    # TODO: REMOVE TOO-LIGHT COLOURS
    col = colmap(float(num)/13.0)
    return col

def Linestyle(simname):
    run, num = RunNum(simname)
    lines = {0:"-",1:"--",2:":"}
    return lines[(num-1) % 3]

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
    print("Oops, something went wrong in linestyles.Labels for simname",simname)
    raise ValueError

#def Legend(run,scatter=True):
#    if run == "imf":
        
