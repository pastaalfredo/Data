#!/usr/bin/env python3

import json, math

# Solve the Poisson equation
def plot_rho(filenm:str, nelec:int, x:list, y:list):
    with open(filenm, "w") as outf:
        dtot = 0
        for i in range(3, len(x)-3):
            d0 = (y[i]-y[i-2])/(x[i]-x[i-2])
            d1 = (y[i+2]-y[i])/(x[i+2]-x[i])
            dd = (d1-d0)/(x[i+1]-x[i-1])
            rho  = -dd/(4*math.pi)
            outf.write("%10g  %10g\n" % ( x[i], rho/nelec ))
            dtot += rho*(x[i]-x[i-1])
        print("Total charge for %s %g" % ( filenm, dtot ))

if __name__ == "__main__":
    T = 100
    elec = { "Br": 36, "F": 10, "Cl": 18, "Li": 2, "Na": 10, "K": 18 }
    with open(f'output_4_{T}.json', 'r') as json_f:
        data = json.load(json_f)
        for elem in data["data"].keys():
            filenm = "rho-" + elem + ".xvg"
            plot_rho(filenm, 1, #elec[elem],
                     data["data"][elem][0], data["data"][elem][1])

    
