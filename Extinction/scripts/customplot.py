'''
Customise plots by importing this
Sam Geen, October 2014
'''

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
plt.ioff()

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

mpl.rc('font', **font)
mpl.rcParams['lines.linewidth'] = 3.0
mpl.rcParams['lines.solid_capstyle'] = "butt"
#mpl.rcParams["savefig.bbox"] = "tight"

# Type 3 Fonts Fix
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.unicode'] = True # Deprecated in Python 3

import fixpyplotraster

#mpl.rc('text',usetex=True)
#preamble = '''
#\usepackage{xcolor}
#\usepackage{color}
#\usepackage{pifont}
#'''
#mpl.rc('text.latex', preamble=preamble)
