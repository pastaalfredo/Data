#!/usr/bin/env python3

import os, math, scipy, sys

def thole(r:float, aiw:float)->float:
    return 1 - (1 + r / (2 * aiw)) * math.exp(- r / aiw)

def bisect(r:float, siw:float):
    a = 1e-6
    b = 10.0
    ta = thole(r, a)-siw
    tb = thole(r, b)-siw
    if ta*tb > 0:
        sys.exit("Cannot solve for siw = %g ta = %g tb = %g" % ( siw, ta, tb ))
    toler = 1e-12
    while (b-a > toler):
        mid = (a+b)/2
        tm  = thole(r, mid)-siw
        if ta*tm < 0:
            b  = mid
            tb = tm
        else:
            a  = mid
            ta = tm
    return (a+b)/2
    
dims = { 
    "LiF":  { "r": 1.64, "sapt": -826.447, "fn": "lithium-fluoride" }, 
    "LiCl": { "r": 2.06, "sapt": -648.698, "fn": "lithium-chloride" }, 
    "LiBr": { "r": 2.26, "sapt": -591.752, "fn": "lithium-bromide" },
    "NaF":  { "r": 2.02, "sapt": -708.783, "fn": "sodium-fluoride" },
    "NaCl": { "r": 2.48, "sapt": -573.763, "fn": "sodium-chloride" },
    "NaBr": { "r": 2.54, "sapt": -562.65, "fn": "sodium-bromide" },
    "KF":   { "r": 2.26, "sapt": -667.864, "fn": "potassium-fluoride" },
    "KCl":  { "r": 2.76, "sapt": -540.113, "fn": "potassium-chloride" },
    "KBr":  { "r": 2.82, "sapt": -538.126, "fn": "potassium-bromide" }
}

kkk = 1389
tex = "ions.tex"
with open(tex, "w") as outf:
    outf.write("""\\begin{table}[ht]
\\centering
\\caption{Electrostatic energy for ion pairs close to the minimum energy distance computed using point charges (PC) and for electrostatics as computed using symmetry-adapted perturbation theory at the SAPT2+3(CCD)$\\delta$MP2 level of theory~\\cite{Parker2014a} with the aug-cc-pvtz basis set~\\cite{Dunning2000a}, computed using the Psi4 suite of programs~\\cite{Psi4}. Details for computing $a$, $\\zeta$ and $\\zeta'$ are given in the running text.}
\\begin{tabular}{lcccccc}
\\hline
 & r  & E$_{PC}$ & E$_{SAPT}$ &  $a$ & $\\zeta$ & $\\zeta'$\\\\
         & (nm) & \\multicolumn{2}{c}{(kJ/mol)} & (nm) & \\multicolumn{2}{c}{(1/nm)} \\\\
\\hline
""")
    for dim in dims:
        vpc   = -kkk/dims[dim]["r"]
        ratio = dims[dim]["sapt"]/vpc
        beta  = "-"
        a     = "-"
        betap = "-"
        rr    = 0.1*dims[dim]["r"]
        if ratio < 1:
            beta  = ( "%.2f" % ( scipy.special.erfinv(ratio)/rr ) )
            aaa   = bisect(rr, ratio) 
            a     = ( "%.5f" % aaa )
            betap = ( "%.2f" % ( 2/(3 * aaa * math.sqrt(math.pi)) ) )
        outf.write("%s & %.4f & %.1f & %.1f & %s & %s & %s\\\\\n" %
                   ( dim, rr, vpc,
                     dims[dim]["sapt"],
                     a, beta, betap ) )

    outf.write("""\\hline
\\end{tabular}
\\label{tab:ionpairs}
\\end{table}
""")

print("Please check %s" % tex)


    
