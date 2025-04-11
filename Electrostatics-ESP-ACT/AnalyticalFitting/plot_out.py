#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import json
from enum import Enum
from potential_elec_functions import *

atnum = { "Li": 3, "F": 9, "Na": 11, "Cl":17, "K": 19, "Br": 35 }

def plot(T:int):
    with open(f'output_4_{T}.json', 'r') as json_f:
        output_data = json.load(json_f)
    data   = output_data['data']

    fig, (axes1, axes2, axes3, axes4, axes5, axes6) = plt.subplots(len(data), 1, figsize=(6, 12))
    axes = [axes1, axes2, axes3, axes4, axes5, axes6]
    clist = [ "F", "Cl", "Br", "Li", "Na", "K" ]
    for i in range(len(clist)):
        compound = clist[i]
        distance_data  = data[compound][0]
        potential_data = data[compound][1]
        if compound=='F' or compound=='Cl' or compound=='Br' or compound=='I':
            charge=-1
        elif compound=='Li' or compound=='Na' or compound=='K':
            charge=1
        else:
            sys.exit("Unknown compound %s" % compound)
        label = f' {compound}'
        if charge == -1:
            label += "-"
        else:
            label += "+"
        for d in range(len(distance_data)):
            if distance_data[d] > 0:
                potential_data[d] -= one_4pi_eps0*charge/distance_data[d]
        axes[i].plot(np.array(distance_data), np.array(potential_data), label=None)
        axes[i].text(.82, .89, label, transform=axes[i].transAxes,  va='top', fontsize=18)
    axes[3].set_ylabel('Residual electrostatic potential (kJ/mol e)', fontsize=18)
    axes[5].set_xlabel(f'Distance ($\AA$)', fontsize=18)
    axes[5].tick_params(axis='x', labelsize=14)
    for ix in range(6):
        axes[ix].tick_params(axis='y', labelsize=14)
    axes1.legend(bbox_to_anchor=(.05, .95), loc='upper left',fontsize=14)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)

    plt.savefig('Test-%d.pdf' % T)
    plt.show()

def plotall(T:int):
    with open(f'output_4_{T}.json', 'r') as json_f:
        output_data = json.load(json_f)
    data   = output_data['data']

    fig, axes = plt.subplots(1, 1)
    clist = [ "F", "Cl", "Br", "Li", "Na", "K" ]
    for i in range(len(clist)):
        compound = clist[i]
        distance_data  = data[compound][0]
        potential_data = data[compound][1]
        if compound=='F' or compound=='Cl' or compound=='Br' or compound=='I':
            charge=-1
        elif compound=='Li' or compound=='Na' or compound=='K':
            charge=1
        else:
            sys.exit("Unknown compound %s" % compound)
        mylabel = f' {compound}'
        if charge == -1:
            mylabel += "-"
        else:
            mylabel += "+"
        for d in range(len(distance_data)):
            if distance_data[d] > 0:
                potential_data[d] -= one_4pi_eps0*charge/distance_data[d]
        axes.plot(np.array(distance_data), np.array(potential_data), label=mylabel)
#        axes.text(.82, .89, label, transform=axes.transAxes,  va='top', fontsize=18)
    axes.set_ylabel('Residual ESP (kJ/mol e)', fontsize=18)
    axes.set_xlabel(f'Distance ($\AA$)', fontsize=18)
    axes.tick_params(axis='x', labelsize=14)
    axes.tick_params(axis='y', labelsize=14)
    axes.legend(loc='upper right',fontsize=18)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)

    plt.savefig('Test-%d.pdf' % T)
    plt.show()

if __name__ == "__main__":
    for T in [ 10, 100 ]:
        plotall(T)
    
