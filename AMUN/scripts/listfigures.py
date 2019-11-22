'''
This gives a list of figures actually used by the paper
Helps speed up run times of scripts by ignoring non papers
Sam Geen, November 2019
'''

def makelist():
    try:
        f = open("figurelist.txt","r")
    except:
        return []
    lines = f.readlines()
    fignames = []
    for line in lines:
        figname = line[line.rfind("/")+1:]
        figname = figname.replace("\n","")
        fignames.append(figname)
    return fignames

if __name__=="__main__":
    print makelist()
