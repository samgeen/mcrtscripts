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
mpl.rcParams['lines.linewidth'] = 2.0
mpl.rcParams['lines.solid_capstyle'] = "butt" # hehehe
mpl.rcParams["savefig.bbox"] = "tight"

# FUCKING TYPE 3 FONTS FUCK FUCK FUCK
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True 

import fixpyplotraster

#mpl.rc('text',usetex=True)
#preamble = '''
#\usepackage{xcolor}
#\usepackage{color}
#\usepackage{pifont}
#'''
#mpl.rc('text.latex', preamble=preamble)
