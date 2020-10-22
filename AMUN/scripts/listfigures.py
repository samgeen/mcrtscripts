'''
This gives a list of figures actually used by the paper
Helps speed up run times of scripts by ignoring non papers
Sam Geen, November 2019
'''

filenames = ["figurelist.txt"]

def makelist(wholepath=False):
    fignames = []
    for filename in filenames:
        try:
            f = open(filename,"r")
        except:
            return []
        lines = f.readlines()
        for line in lines:
            figname = line
            if not wholepath:
                figname = figname[figname.rfind("/")+1:]
            figname = figname.replace("\n","")
            fignames.append(figname)
    return fignames

if __name__=="__main__":
    print(makelist())
