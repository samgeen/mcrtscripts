'''
Make a table of all the R values by plot
'''

from startup import *
import glob

def rstr(rin):
    fr = float(rin)
    rout =  "%.3g" % fr
    if len(rout[rout.find("."):]) == 3:
        rout = rout+r"0"
    if rout[0] != "-":
        #rout = r"\hphantom{-}"+rout
        rout = r"+"+rout
    return rout

def run(relation,params,labels):
    nl = os.linesep
    print "Making table for", relation
    intro = r'''
    \begin{center}
    \begin{tabular}{llllllll}
    \textbf{Quantity} & \multicolumn{3}{|c|}{\textbf{Correlation coefficient}} \\
    '''
    # Comparison didn't work super well for IC runs here so I left it out
    #if relation == "starrelations":
    #    intro += r"Run = & \textsc{stars} & \textsc{turb} \\ "+os.linesep+" "
    if relation == "correlatestructure":
        intro += r"$n_{\mathrm{H}}\mathrm{(cutoff)/cm}^{-3} = $ & 10 & 100 & 1000 \\ "+os.linesep+" "
    if relation == "correlatestructure_time":
        intro += r"$t / t_{ff} = $  & 0.5 & 1 & 1.5 & 2 \\ "+os.linesep+" "
    intro += r"\hline "+os.linesep+" "
    outro = r'''
    \end{tabular}
    \end{center}
    '''

    table = intro+""
    for param,label in zip(params,labels):
        if relation == "starrelations":
            fnameic = "../plots/"+relation+"_"+param+"ic_rxy.txt"
            fnameimf = "../plots/"+relation+"_"+param+"imf_rxy.txt"
            rvalueic = rstr(open(fnameic,"r").read().strip())
            rvalueimf = rstr(open(fnameimf,"r").read().strip())
            #table += label+" & "+rvalueimf+" & "+rvalueic+" \\\\ "+os.linesep+" "
            table += label+" & "+rvalueimf+" \\\\ "+os.linesep+" "
        elif relation == "correlatestructure":
            table += label+" & "
            first = True
            for n in ["10p0","100p0","1000p0"]:
                if first:
                    first = False
                else:
                    table += " & "
                fname = "../plots/"+relation+"_"+param+"_n"+n+"parameter_rxy.txt"
                rvalue = rstr(open(fname,"r").read().strip())
                table +=rvalue+" "
            table += "\\\\ "+os.linesep+" "
        else:
            table += label+" & "
            first = True
            n = "10p0"
            for t in ["0p5","1p0","1p5","2p0"]:
                if first:
                    first = False
                else:
                    table += " & "
                fname = "../plots/correlatestructure_"+param+"_n"+n+"_"+t+"tff_parameter_rxy.txt"
                rvalue = rstr(open(fname,"r").read().strip())
                table +=rvalue+" "
            table += "\\\\ "+os.linesep+" "
    table += outro
    print "~~~"
    print table
    print "~~~"
    f = open("../plots/table_"+relation+".tex","w")
    f.write(table)
    f.close()

if __name__=="__main__":
    run("correlatestructure",
        ["T","L","M","S","C"],
        ["Triaxiality","Longest Axis","Middlemost Axis","Shortest Axis","Compactness"])
    run("starrelations",
        ["alltimemax","compactness","firstmass","nphotons","nphotonstot","nphotonstff","tracklength"],
        ["Most Massive Star","Cluster Compactness","Mass of First Star Formed","Peak Photon Emission Rate",
         "Total Photons Emitted By Cluster","Photons Emitted in One Freefall Time","Distance Travelled By Star"])
    run("correlatestructure_time",
        ["T","L","M","S","C"],
        ["Triaxiality","Longest Axis","Middlemost Axis","Shortest Axis","Compactness"])

