#!/usr/bin/env python3
import os

npoints = 10
start_r = 1.8
step = 0.03

for i in range(npoints):
    r = start_r + i * step
    tag = f"{i:04d}"

    if not os.path.exists(tag):
        os.makedirs(tag)

    py_filename = os.path.join(tag, f"{tag}.py")
    slurm_filename = os.path.join(tag, f"{tag}.slurm")

    with open(py_filename, "w") as f:
        f.write(f"""#!/usr/bin/env python3
#SBATCH -t 32:00:00
#SBATCH -c 4
#SBATCH -p CLUSTER,CLUSTER-AMD
import os, sys
import numpy as np
import psi4 as psi4
psi4.core.set_num_threads(4)
psi4.set_options({{
    'guess': 'read',
    'reference': 'rhf',
    'cachelevel': 1,
    'print': 1,
    'damping_percentage': 20
}})
psi4.set_memory(12000000000)
psi4_io = psi4.core.IOManager.shared_object()
tmpdir = "/scratch"
psi4_io.set_default_path(tmpdir)
psi4.core.set_output_file('{tag}.log', False)

geometry = \"\"\"
 1 1
 Li 0.0 0.0 0.0
 --
 -1 1
 Cl 0.0 0.0 {r:.16f}
\"\"\"
geom = psi4.geometry(geometry)
psi4.basis_helper(\"\"\"
assign aug-cc-pvtz              
\"\"\")

forces = None
mydict = {{}}
mydict["energies"] = {{}}
mydict["energies"]["energy"] = psi4.energy("sapt2+3(ccd)dmp2")
mydict["energies"]["Electrostatics"] = psi4.variable('SAPT ELST ENERGY')
mydict["energies"]["Exchange"] = psi4.variable('SAPT EXCH ENERGY')
mydict["energies"]["Induction"] = psi4.variable('SAPT IND ENERGY')
mydict["energies"]["Dispersion"] = psi4.variable('SAPT DISP ENERGY')
mydict["energies"]["InteractionEnergy"] = psi4.variable('SAPT TOTAL ENERGY')

with open("{tag}.out", "w") as result:
    result.write("    1\\n")
    result.write(" Li 0.0 0.0 0.0\\n")
    result.write("    1\\n")
    result.write(" Cl 0.0 0.0 {r:.16f}\\n")
    result.write("Energy %12.8f \\n" % mydict["energies"]["energy"])
    if None != forces and None != forces.all():
        for f in forces:
            result.write("%12.8f  %12.8f  %12.8f\\n" % (f[0], f[1], f[2]))
""")

    with open(slurm_filename, "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -J sapt_{tag}
#SBATCH -o {tag}.slurm.out
#SBATCH -e {tag}.slurm.err
#SBATCH -t 32:00:00
#SBATCH -c 4
#SBATCH -p CLUSTER,CLUSTER-AMD
module load psi4
python3 {tag}.py
""")

    os.system(f"cd {tag} && sbatch {tag}.slurm")
