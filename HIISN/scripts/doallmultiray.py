'''
Run all multirays
Sam Geen, January 2016
'''

import Hamu
import multiray
import linestyles # contains sim names

def findtsn(ws):
    if ws == "HIISN":
        return 5.52990445
    elif ws == "HIISN4":
        return 4.25

workspaces = ["HIISN","HIISN4"]

for ws in workspaces:
    Hamu.Workspace(ws)
    tsn = findtsn(ws)
    tff = tsn-3.0
    ems = linestyles.makeems()
    sns = ["-NSN","-SN"]
    # Compare each run at a given emission rate
    for em in ems:
        if not "N50" in sns:
            simnames = [em+sn for sn in sns]
        else:
            simnames = [em+sn for sn in ["-NSN","-SN","-HN","-MSN"]]
        times = [tff,tff+1.0,tff+2.0,tsn,tsn+0.5,tsn+1.0]
        multiray.PlotForSims(simnames,times,name=em)
