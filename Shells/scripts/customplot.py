'''
Customise plots by importing this
Sam Geen, October 2014
'''

import matplotlib as mpl
mpl.use("Agg")

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

mpl.rc('font', **font)
mpl.rcParams['lines.linewidth'] = 3.0
mpl.rcParams['lines.solid_capstyle'] = "butt" # I have no idea why it's called this
mpl.rcParams["savefig.bbox"] = "tight"

# Fixes some problems I had submitting with "Type 3 fonts", whatever those are
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True 

import fixpyplotraster