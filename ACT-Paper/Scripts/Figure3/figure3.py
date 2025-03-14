#!/usr/bin/env python3

import os

myffs = { "OPLS": "alexandria train_ff -nooptimize -ff ../../../ForceField/OPLS2020-water-alcohol.xml -charges ../../../MolProps/OPLS2020-charges.xml -mp ../../../MolProps/sapt-alcohol-HH-OPLS.xml -sel ../../../Selection/alcwater.dat",
          "model-BC": "alexandria train_ff -nooptimize -ff ../../../ForceField/model-BC.xml -charges ../../../MolProps/mp2-alcohol-HH.xml -mp ../../../MolProps/sapt-alcohol-HH.xml -sel ../../../Selection/alcwater.dat"
         }

for ff in myffs:
    os.makedirs(ff, exist_ok=True)
    os.chdir(ff)
    os.system(myffs[ff])
    os.chdir("..")

for xvg in [ "EPOT", "COULOMB", "EXCHANGE", "DISPERSION" ]:
    pdf = ( "%s.pdf" % xvg )
    cmd = ( "viewxvg -sq -ls None -mk -pdf %s -noshow -alfs 24 -lfs 24 -tickfs 24 -f" % ( pdf ) )
    for ff in myffs:
        cmd += ( " %s/%s.xvg " % ( ff, xvg ) )
    cmd += " -label"
    for ff in myffs:
        cmd += ( " %s " % ( ff ) )
    os.system(cmd)
    print("Please check %s" % pdf)

    
