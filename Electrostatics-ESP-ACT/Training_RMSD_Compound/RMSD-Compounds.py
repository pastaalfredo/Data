#!/usr/bin/env python3

import os

def run_one(filenm: str) -> dict:
    tempf = "rmsd.txt"
    os.system("grep 'COULOMB RMSD' %s > %s" % (filenm, tempf))
    
    mydict = {}
    with open(tempf, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if len(words) == 14:
                mydict[words[13]] = {"RMSD": words[2], "N": words[12]}
    return mydict


def run_legacy(filenm: str) -> dict:
    tempf = "rmsd_legacy.txt"
    os.system("grep 'COULOMB RMSD' %s > %s" % (filenm, tempf))
    
    mydict = {}
    with open(tempf, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if len(words) == 14:
                mydict[words[13]] = words[2]
    return mydict


if __name__ == "__main__":
    basefn = "train-Inter6.log"
    alldata = {}
    for mydir in ["PC", "GC", "PC+GV", "PC+GVS"]:
        alldata[mydir] = run_one("Elec/" + mydir + "/" + basefn)
    
    legacy_files = {
        "qMulliken2": "legacy/Mulliken.log",
        "qHirshfeld2": "legacy/Hirshfeld.log",
         "qESP2": "legacy/ESP.log",
         "qCM52": "legacy/CM5.log"
    }
    legacy_data = {key: run_legacy(path) for key, path in legacy_files.items()}
    
    texfn = "rmsdtable.tex"
    with open(texfn, "w") as outf:
        outf.write("\\begin{longtable}{lcccccccccc}\n")
        outf.write("\\caption{Root mean square deviation (kJ/mol) from SAPT2+(CCD)$\delta$MP2 electrostatics per compound dimer for the different ACT models and widely-used models, Mulliken, Hirshfeld, ESP, and CM5. N is the number of conformations of each dimer used.}\\\\\n")
        outf.write("\\hline\n")
        outf.write("Dimer & N ")
        for md in alldata.keys():
            outf.write(" & %s" % md)
        for legacy in legacy_data.keys():
            outf.write(" & %s" % legacy)
        outf.write("\\\\\n")
        outf.write("\\hline\n")
        for dimer in sorted(alldata["PC"].keys()):
            outf.write("%s & %s " % (dimer.replace("#", "-"), alldata["PC"][dimer]["N"]))
            for legacy in legacy_data.keys():
                outf.write(" & %s" % legacy_data[legacy].get(dimer, "-"))
            for md in alldata.keys():
                outf.write(" & %s" % alldata[md][dimer]["RMSD"])
            outf.write("\\\\\n")
        outf.write("\\hline\n")
        outf.write("\\end{longtable}\n")
