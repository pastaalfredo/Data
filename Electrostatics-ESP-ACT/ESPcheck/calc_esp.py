#!/usr/bin/env python3

import os

def run_one(ff:str)->float:
    myff = os.path.basename(ff)
    mylog = myff.replace(".xml", ".log")
    os.system("alexandria train_ff -ff ../AlexandriaFF/%s -mp ../../../ACTdata/MolProps/hf-aug-cc-pvtz.xml -nooptimize -sel ../Selection/monomer.dat -g %s -v 4" % ( myff, mylog ) )
    rmsd = None
    with open(mylog, "r") as inf:
        for line in inf:
            if line.find("ESP (kJ") >= 0:
                words = line.split()
                try:
                    rmsd = float(words[6])
                except ValueError:
                    print("Cannot read ESP RMSD in %s" % mylog)
    return rmsd
    
if __name__ == "__main__":
    allffs = [ "all-p.xml", "coul-gv.xml", "esp-gv.xml", "all-g.xml",
    	       "all-pg.xml", "coul-p.xml", 
               "all-gv.xml", "coul-g.xml", "esp-g.xml" ]
    myrmsd = {}
    for ff in allffs:
        print("Will run %s" % ff)
        rmsd = run_one(ff)
        print("==================")
        if rmsd:
            myrmsd[ff[:-4]] = rmsd
    for key in sorted(myrmsd.keys()):
        val = myrmsd[key]
        print("%s  %g" % ( key, val ))

