'''
Make input tables for cooling stuff for Eric
Sam Geen, January 2021
'''

from startup import *

import pickle as pik

def maketable(snap,hydros):
    '''
    Make an ascii table of every cell value
    '''
    amr = hydrofuncs.amr_source(snap,hydros) 
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    vars = {}
    for hydro in hydros:
        vars[hydro] = hydrofuncs.scale_by_units(snap, hydro)(cells)
    return vars

def MakeLabels(hydros):
    labels = {}
    for hydro in hydros:
        labels[hydro] = hydrofuncs.hydro_label(hydro)
    return labels

def runforsim(simname,hydros):
    sim = Hamu.Simulation(simname)
    for snap in sim.Snapshots():
        vars = maketable(snap,hydros)
        labels = MakeLabels(hydros)
        output = str(snap.OutputNumber().zfill(5)
        f = open("../coolinginputs/"+simname+"_"+output+".pik","wb")
        pik.dump(vars,f)
        pik.dump(labels,f)
        f.close()

def run():
    hydros = ["dx","T","nH","xHIII","xHeII","xHeIII","NpHIII","NpHeII","NpHeIII"]
    runforsim("UVWINDPRESS_30",hydros)

if __name__=="__main__":
    run()