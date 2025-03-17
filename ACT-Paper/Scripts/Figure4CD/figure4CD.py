#!/usr/bin/env python3

import os, json

myffs = { "OPLS": "alexandria nma -ff ../../../FF_ACT/OPLS2020-water-alcohol.xml -charges ../../../MolProps/OPLS2020-charges.xml -json ",
          "model-C": "alexandria nma -ff ../../../FF_ACT/model-C.xml -charges ../../../MolProps/mp2-alcohol-HH.xml -json "
         }

for ff in myffs:
    os.makedirs(ff, exist_ok=True)
    os.chdir(ff)
    os.system(myffs[ff] + (" -db ethanol -ir IRethanol ") )
    os.chdir("..")


xmin = 0
xmax = 3800
mol = "ethanol"
pdf = "figure4C.pdf"
os.system("viewxvg -f OPLS/IR%s.xvg model-C/IR%s.xvg -label 'OPLS2020' 'Model C' -lfs 28 -alfs 28 -tickfs 24 -pdf %s -noshow  -xmin %d -xmax %d -colors orange green" % ( mol, mol, pdf, xmin, xmax ) )
print("Please check %s" % pdf)
    
Exper = "Experiment"
freqs = { Exper : [ 200, 251, 417, 812, 888, 1028, 1091, 1161,
                    1256, 1275, 1371, 1412, 1446, 1464, 1490,
                    2900, 2910, 2939, 2984, 2991, 3653 ] }
for ff in myffs:
    os.chdir(ff)
    with open("nma.json", "r") as inf:
        mydata = json.load(inf)
        freqs[ff] = []
        index = 1
        for fq in mydata["simulate"][1]["Frequencies"][0]["Alexandria"]:
            freqs[ff].append(fq[str(index)][0]["value"])
            index += 1
    os.chdir("..")

freqxvg = "allfreqs.xvg"
with open(freqxvg, "w") as outf:
    outf.write("@ xaxis label 'Experimental Frequency (1/cm)'\n")
    for f in range(len(freqs[Exper])):
        outf.write("%10g" % freqs[Exper][f])
        for ff in myffs:
            outf.write("  %10s" % freqs[ff][f])
        outf.write("\n")
pdf = "figure4D.pdf"
os.system("viewxvg -f %s -pdf %s -noshow -res -ls None -mk -label OPLS 'Model C' -lfs 28 -alfs 28 -tickfs 24" % ( freqxvg, pdf ) )
print("Please check %s" % pdf)
