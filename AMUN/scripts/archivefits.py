'''
Archive the RDM files used in the paper
Sam Geen, September 2020
'''

from startup import *

import os, sys

import listfigures

def runfits():
    '''
    Copy fits files from plotting to an archive file
    '''
    fignames = listfigures.makelist(wholepath=True)
    for figname in fignames:
        fitsname = figname.replace(".pdf",".fits")
        os.system("cp "+fitsname+" ../archive/")
    os.system("zip -r ../plots/archive.zip ../archive/")

def makedir(folder):
    try:
        os.mkdir(folder)
    except:
        pass
    
def runsims():
    '''
    Copy simulation files
    '''
    masses = ["30","60","120","120_DENSE"]
    simnames = ["NOFB","NOFB_DENSE"]
    simnames += ["UV_"+mass for mass in masses]
    simnames += ["UVWIND_"+mass for mass in masses]
    simnames += ["UVWINDPRESS_"+mass for mass in masses]
    for simname in simnames:
        simfolder = simfolders[simname]
        rdmfolder = "../RDM/"+simname
        makedir(rdmfolder)
        copylist = ["*3d","job*sh","*.nml","pymses_field_descrs.py"]
        for tocopy in copylist:
            os.system("cp "+simfolder+"/"+tocopy+" "+rdmfolder)
        # Copy ramses.data file separately
    os.system("cp "+simfolder+"/ramses.data ../RDM")
    
if __name__=="__main__":
    runsims()
