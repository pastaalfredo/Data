%mem=12000MB
%nprocshared=8
%chk=methylammonium19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB

methylammonium

1 1
N 0.000000 -0.000000 0.712271
C -0.000000 0.000000 -0.803230
H 0.000000 0.953575 1.088129
H 0.825820 -0.476787 1.088129
H -0.825820 -0.476787 1.088129
H 0.894550 0.516469 -1.143633
H -0.894550 0.516469 -1.143633
H -0.000000 -1.032937 -1.143633


