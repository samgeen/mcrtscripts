import findproperties
from pymses.filters import CellsToPoints

def readfunc(funcname,lines):
    newfunc = ""
    copyline = False
    indent = 0
    linecount = 0
    for line in lines:
        #print line
        # Check for correct function definition
        if "def "+funcname in line:
            copyline = True
            #indent = line.find("def ")
        # Check for indent size
        if linecount == 1:
            indent = line.find(line.strip())
        # Check for an empty line
        if len(line.strip()) < indent:
            copyline = False
        if copyline:
            linecount += 1
            newfunc += line+"\n"
    # Add warning we've modified the function
    newlines = newfunc.split("\n")
    newlines.insert(1," "*indent+'print "WARNING: CODE MODIFIED BY skimfunc.py"')
    newfunc = "\n".join(newlines)
    return newfunc.strip()

def findamrs(funccode):
    newfunc = ""
    lines = funccode.split("\n")
    search = "amr_source("
    search2 = "cells = CellsToPoints(amr).flatten()"
    amrs = []
    for line in lines:
        if search in line:
            # Get list of hydro variables we want
            newamrs = line[line.find(search)+len(search):]
            newamrs = newamrs[:newamrs.find(")")]
            amrs += eval(newamrs)
            # Replace amr_source call with full list pre-computed
            newamrs = line[line.find(search):]
            line = line.replace(newamrs,"full_amr")
        if search2 in line:
            line = line.replace(search2,"cells = snap.full_cells")
        newfunc += line+"\n"
    return amrs, newfunc.strip()
    
def skimfuncs(funcs, filename):
    '''
    Skim functions to make new code to run
    '''
    funccode = {}
    amrs = []
    print "Skimming functions for pre-computing...",
    # Open file
    lines = open(filename).read().split("\n")
    # Make new functions that use the full_amr function
    for func in funcs:
        funccode[func] = readfunc(func,lines)
        newamrs, newfunc = findamrs(funccode[func])
        funccode[func] = newfunc
        amrs += newamrs
    # Find unique elements
    amrs = list(set(amrs))
    print "Done!"
    return amrs, funccode

def makefullamr(snap,amrs):
    '''
    Make a full AMR source object and attach it to the snapshot
    '''
    full_amr = snap.amr_source(amrs)
    snap.full_amr = full_amr
    snap.full_cells = CellsToPoints(full_amr).flatten()
    return snap

# Fake snapshot for testing
class Snapshot(object):
    def __init__(self):
        self._a = 4.0

    def amr_source(self,amrs):
        print "Running amr_source dummy function"
        return "FULL AMR"


if __name__=="__main__":
    funcs = ["testfunc","secondfunc"]
    filename = "findproperties.py"
    amrs, funccode = skimfuncs(funcs,filename)
    print amrs
    snap = Snapshot()
    makefullamr(snap,amrs)
    for func, code in funccode.iteritems():
        print func+":"
        #print code
        exec(code)
    testfunc(snap)
    secondfunc(snap)
