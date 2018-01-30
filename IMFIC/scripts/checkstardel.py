'''
Check for a star being deleted
Sam Geen, January 2018
'''

from startup import *

f = open("../plots/deletelog.txt","w")

def runforsim(simname):
    sim = hamusims[simname]
    # Get first deletion time
    folder = sim.Folder()
    log = folder+"/LOG"
    logtxt = open(log,"r").read()
    firstdel = logtxt[:logtxt.find("DELETE")]
    time = firstdel[firstdel.rfind(" t= "):]
    time = time[:time.find("dt")]
    t = float(time[time.find("t=")+2:])
    snap = sim.FindAtTime(t,forceBefore=True)
    lastout = snap.OutputNumber()
    nl = "\n"
    f.write( simname +nl)
    f.write("First stellar deletion at "+str(t)+nl)
    f.write( "Last snapshot before first stellar deletion: "+str(lastout) +nl )
    f.write("---"+nl)
    # Purge later outputs
    for snap in sim.Snapshots():
        if snap.OutputNumber() > lastout:
            trash = folder+"/trash/"
            try:
                os.makedirs(trash)
            except:
                pass
            cmd = "mv "+folder+"/output_"+str(snap.OutputNumber()).zfill(5)+" "+trash
            print cmd
            os.system(cmd)

def run(simnames,plotname):
    for simname in simnames:
        runforsim(simname)

if __name__=="__main__":
    run(imfsims,"imf")
    run(icsims,"ic")
    f.close()
