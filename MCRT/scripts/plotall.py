'''
Plot all quickplots for a suite
Sam Geen, July 2014
'''

import quickplots

suites = {"nomhd": "03_diffstars_lvl10",
          "mhd": "04_diffstars_mhd10"}
sims = {"NoRT":"00_nort",
            "1e49":"01_1e49_O5V",
            "1e48":"02_1e48_B0V",
            "1e47":"03_1e47_B1V"}

if __name__=="__main__":
    for sname, sfolder in suites.iteritems():
        for name, folder in sims.iteritems():
            quickplots.Run("../runs/"+sfolder+"/"+folder,
                           name+"_"+sname)
