'''
Make a table of all the R values by plot
'''

from startup import *
import glob

def run(relation,params,labels,caption):
    nl = os.linesep
    print "Making table for", relation
    intro = r'''
    \begin{table*}
    \begin{center}
    \begin{tabular}{llllllll}
    \textbf{Quantity} & \textbf{Correlation coefficient} \\
    '''
    if relation == "correlatestructure":
        intro += r"$n_{\mathrm{H}}\mathrm{(cutoff)/cm}^{-3} = $ & 10 & 100 & 1000 \\ "+os.linesep+" "
    intro += r"\hline "+os.linesep+" "
    outro = r'''
    \end{tabular}
    \end{center}
    \caption{CAPTION}
    \label{Rtable} 
    \end{table*}
    '''.replace("CAPTION",caption)

    table = intro+""
    for param,label in zip(params,labels):
        if relation == "starrelations":
            fname = "../plots/"+relation+"_"+param+"_rxy.txt"
            rvalue = open(fname,"r").read().strip()
            table += label+" & "+rvalue+" \\\\ "+os.linesep+" "
        else:
            table += label+" & "
            first = True
            for n in ["10p0","100p0","1000p0"]:
                if first:
                    first = False
                else:
                    table += " & "
                fname = "../plots/"+relation+"_"+param+"_n"+n+"parameter_rxy.txt"
                rvalue = open(fname,"r").read().strip()
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
        ["Triaxiality","Longest Axis","Middlemost Axis","Shortest Axis","Compactness"],
        "CAPTION HERE")
    run("starrelations",
        ["alltimemax","compactness","firstmass","nphotons","nphotonstot"],
        ["Most Massive Star","Cluster Compactness","Mass of First Star Formed","Peak Photon Emission Rate",
         "Total Photons Emitted By Cluster"],
        "CAPTION GOES HERE")
