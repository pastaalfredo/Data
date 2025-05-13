#!/usr/bin/env python3
#SBATCH -t 32:00:00
#SBATCH -c 4
#SBATCH -p CLUSTER,CLUSTER-AMD
import os, sys
import numpy as np
import psi4 as psi4
psi4.core.set_num_threads(4)
psi4.set_options({'guess': 'read', 'reference': 'rhf', 'cachelevel': 1, 'print': 1, 'damping_percentage': 20})
psi4.set_memory(12000000000)
psi4_io = psi4.core.IOManager.shared_object()
tmpdir = "/scratch"
psi4_io.set_default_path(tmpdir)
psi4.core.set_output_file('0259.log', False)
geometry= """
 1 1
 Li 0.0 0.0 0.0
 --
 -1 1
 Cl 0.0 0.0 9.31999999999989
"""
geom = psi4.geometry(geometry)
psi4.basis_helper("""
assign aug-cc-pvtz
""")
forces = None
mydict = {}
mydict["energies"] = {}
mydict["energies"]["energy"] = psi4.energy("sapt2+3(ccd)dmp2")
mydict["energies"]["Electrostatics"] = psi4.variable('SAPT ELST ENERGY')
mydict["energies"]["Exchange"] = psi4.variable('SAPT EXCH ENERGY')
mydict["energies"]["Induction"] = psi4.variable('SAPT IND ENERGY')
mydict["energies"]["Dispersion"] = psi4.variable('SAPT DISP ENERGY')
mydict["energies"]["InteractionEnergy"] = psi4.variable('SAPT TOTAL ENERGY')
with open("0259.out", "w") as result:
    result.write("    1\n")
    result.write(" Li 0.0 0.0 0.0\n")
    result.write("    1\n")
    result.write(" Cl 0.0 0.0 9.31999999999989\n")
    result.write("Energy %12.8f \n" % mydict["energies"]["energy"])
    if None != forces and None != forces.all():
        for f in forces:
            result.write("%12.8f  %12.8f  %12.8f\n" % (f[0], f[1], f[2]))
