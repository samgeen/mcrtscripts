'''
Code to run on script startup for this project
Includes imports that can be used by other scripts
Sam Geen, June 2016
'''

import sys, os, glob
import numpy as np
# Import from the main scripts folder
sys.path.append("../../scripts")

import customplot
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import hydrofuncs
import snaptime
import simlabels

import linestyles

import Hamu
import pymses

import numpy as np

Hamu.Workspace("YuleLads")

allsims = simlabels.allsims  

# Physical conversions
X = 0.76
kB = 1.38062e-16
pcincm = 3.086e18
Msuning = 1.9891e33
mHing = 1.66e-24
Myrins = 3.1556926e13

# Some global definitions
Msolar = "M$_{\odot}$"

def MakeDirs(folder):
    try:
        print "Making directory", folder
        os.makedirs(folder)
    except:
        print "Not making directory", folder, "- it already exists!"
        pass

print "Imported various project-global modules"
