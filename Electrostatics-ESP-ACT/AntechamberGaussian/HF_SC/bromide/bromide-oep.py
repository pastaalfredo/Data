#!/usr/bin/env python3

#SBATCH -A naiss2023-5-531
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 72:00:00

import sys, os, openbabel
from gaussian     import *
from actutils     import *
from elements     import *
from get_mol_dict import *
from mol_csv_api  import *
from dbutils      import *

def write_com(M, molname, ilot, atoms, ncores, mem, cont=False):
    mol = M.find_mol(molname)
    if not mol:
        sys.exit("Cannot find %s in alexandria.csv" % molname)
        
    com = ("%s-%d-oep.com" % ( molname, ilot ) )
    with open(com, "w") as outf:
        outf.write("%%mem=%dMB\n" % mem)
        outf.write("%%nprocshared=%d\n" % ncores)
        outf.write("%%chk=%s%d.chk\n" % ( molname, ilot ))
        if ilot == 1:
            outf.write("#P B3LYP/gen Opt=(Redundant, calcall) symm=(loose,follow)\nmaxdisk=128GB\n\n")
        elif ilot == 3:
            if cont:
                outf.write("#P B3LYP/Gen Geom=AllCheckpoint Int(Grid=SuperFineGrid) Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq\nmaxdisk=128GB\n\n")
            else:
                outf.write("#P B3LYP/Gen Int(Grid=SuperFineGrid) Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq\nmaxdisk=128GB\n\n")
        elif ilot == 19:
            outf.write("#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq\nmaxdisk=128GB\n\n")
        
        if not cont:
            outf.write("%s\n\n" % molname)
            outf.write("%d %d\n" % ( mol.charge, mol.mult))
        haveElem = {}
        for i in range(1,1+len(atoms.keys())):
            elem = AtomNumberToAtomName(atoms[i]["atomic_number"])
            if not elem:
                sys.exit("Could not find element for " + str(atoms[i]["atomic_number"]))
            if not cont:
                outf.write("%s %f %f %f\n" % ( elem, atoms[i]["X"], atoms[i]["Y"], atoms[i]["Z"]))
            haveElem[elem] = 1
        if not cont:
            outf.write("\n")
        
        liquids = os.environ["LIQUIDS"]
        specialElements = "I" in haveElem or "Kr" in haveElem or "Ar" in haveElem
        basis = None
        if ilot == 1:
            basis = ( "%s/MOLECULES/BASIS/6-311Gxx" % liquids)
        elif ilot == 3:
            if specialElements:
                basis = ( "%s/MOLECULES/BASIS/augccpVTZ-I" % liquids)
            else:
                basis = ( "%s/MOLECULES/BASIS/augccpVTZ" % liquids)
        if basis:
            with open(basis, "r") as inf:
                for line in inf:
                    outf.write(line)
            outf.write("\n")
        if specialElements:
            for i in range(6):
               outf.write("I, 2.1\n")
               #outf.write("Ar,1.88\n")
               #outf.write("Kr,2.02\n")
               outf.write("\n")
        outf.write("\n")

    return com
    
def run_gauss(com, ncores):
    logfile = com[:-3] + "log"
    status  = run_g16(com, logfile, ncores)
    return status, logfile

def run_both(M, molname, adb, ncores, mem, ilot:int):
    target   = ( "../../../OEP/%s/%s-%d-oep.log.gz" % ( molname, molname, ilot ))
    if os.path.exists(target):
        print("%s already done" % target)
        return
    sdffile  = adb + "/COMPOUNDS/" + molname + ".sdf"
    md       = MoleculeDict()
    md.read(sdffile, "sdf")
    com1     = write_com(M, molname, 1, md.atoms, ncores, mem)
    logfile1 = com1[:-3] + "log"
    status   = Status.unclear
    if os.path.exists(logfile1): 
        status = check_normal_termination(logfile1)
    if not status == Status.normal_termination:
        status, logfile1 = run_gauss(com1, ncores)
    if status == Status.normal_termination:
        md2      = MoleculeDict()
        md2.read(logfile1, "g09")
        com3     = write_com(M, molname, ilot, md2.atoms, ncores, mem)
        logfile3 = com3[:-3] + "log"
        if os.path.exists(logfile3):
            status = check_normal_termination(logfile3)
            if not status == Status.normal_termination:
                if status == Status.unclear:
                    chkfile = com3[:-3] + "chk"
                    cont = os.path.exists(chkfile)
                    com3 = write_com(M, molname, ilot, md2.atoms, ncores, mem, cont)
                    status, logfile3 = run_gauss(com3, ncores)
                else:
                    for fn in [ logfile3, ( "%s.chk" % molname ) ]:
                        if os.path.exists(fn):
                            os.unlink(fn)
                    run_both(M, molname, adb, ncores, mem)
        else:
            status, logfile3 = run_gauss(com3, ncores)
    return status
    
if __name__ == '__main__':
    adb = os.environ["AlexandriaDB"]
    M = Molecules()
    M.read(act_library_filename("alexandria.csv"), 3, False)
    ncores   = 8
    scon     = "SLURM_CPUS_ON_NODE"
    if scon in os.environ:
        ncores   = int(os.environ[scon])
    mem      = 1500*ncores
    ilot     = 19
    status = run_both(M, "bromide", adb, ncores, mem, ilot)
    print(status)
    
