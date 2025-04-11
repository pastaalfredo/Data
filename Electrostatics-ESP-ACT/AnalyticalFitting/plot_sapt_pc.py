#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import json
import math
import os
from potential_elec_functions import *

figs = "Figures"

def main(T:int):
    cations = ["Li", "Na", "K"]
    anions = ["F", "Cl", "Br"]
    dir = 'SAPT_Alkali_Halides'

    fig, allaxes = plt.subplots(3, 3, figsize=(6, 6))
    axes = allaxes.flatten()
    i = 0
    # Start RMSD calc from this distance in Angstrom
    xstart = 2
    istart = 0
    for cation in cations:
        for anion in anions:
            compound = cation + anion
            myfile = f'SAPT_Alkali_Halides/distances_Electrostatics-{cation.lower()}{anion.lower()}.txt'
            dataset = np.loadtxt(myfile)
            x, y = dataset[:, 0], dataset[:, 1]
            sorted_indices = np.argsort(x)
            while (x[sorted_indices][istart] < xstart and
                   istart < len(x[sorted_indices]) - 1):
                istart += 1
            x1 = x[sorted_indices][istart:]
            y1 = y[sorted_indices][istart:]

            axes[i].plot(x1, y1, label=("SAPT %s" % compound), color='r')
            ener_pc = []
            for j in range(len(x1)):
                ener = 0
                if x1[j] > 0:
                    ener = -one_4pi_eps0/x1[j]
                    ener_pc.append(ener)
            # Compute RMSD
            pcrmsd = np.sqrt(np.mean((y1-ener_pc)**2))
            pclabel = ("PC (%.0f kJ/mol)" % ( pcrmsd ))
            #pclabel = "PC"
            axes[i].plot(x1, ener_pc, label=pclabel, color='blue')

            if i == 3:
                axes[i].set_ylabel('Electrostatic energies (kJ/mol)', fontsize=18)
            if i == 7:
                axes[i].set_xlabel('Distance ($\AA$)', fontsize=18)
            axes[i].tick_params(axis='x', labelsize=14)
            axes[i].tick_params(axis='y', labelsize=14)
            axes[i].legend(fontsize=16)
            i += 1
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    plt.show()
#    plt.savefig(f'Figures/SAPT_PC_{T}.pdf')


if __name__ == "__main__":
    os.makedirs(figs, exist_ok=True)
    for T in [ 10 ]:
        main(T)
