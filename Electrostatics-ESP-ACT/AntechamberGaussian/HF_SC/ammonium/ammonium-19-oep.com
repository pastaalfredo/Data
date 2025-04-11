%mem=12000MB
%nprocshared=8
%chk=ammonium19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB

ammonium

1 1
H 0.592548 0.592548 0.592548
N 0.000000 0.000000 0.000000
H -0.592548 -0.592548 0.592548
H -0.592548 0.592548 -0.592548
H 0.592548 -0.592548 -0.592548


