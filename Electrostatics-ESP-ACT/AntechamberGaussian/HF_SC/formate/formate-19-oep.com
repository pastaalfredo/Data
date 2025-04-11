%mem=12000MB
%nprocshared=8
%chk=formate19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB

formate

-1 1
O -1.137576 -0.206870 0.000000
O 1.137597 -0.206755 0.000000
C -0.000016 0.307620 0.000000
H -0.000074 1.463278 0.000000


