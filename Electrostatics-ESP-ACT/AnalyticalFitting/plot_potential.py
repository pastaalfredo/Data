#!/usr/bin/env python3

from potential_elec_functions import *

def plot(outfn, alpha, func):
    with open(outfn, "w") as outf:
        for r in range(0,501):
            x = 0.01*r
            y = eval(func)(x, alpha)
            outf.write("%10f  %10f\n" % ( x, y ))

for alpha in [ 1, 3, 5, 10 ]:
    for func in [ "slater_charge", "slater2_charge", "gaussian" ]:
        outfn = str(func) + "-alpha-" + str(alpha) + ".xvg"
        plot(outfn, alpha, func)

           
    
