%mem=12000MB
%nprocshared=8
%chk=water19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB

water

0 1
H 0.000000 0.756998 -0.474863
O 0.000000 0.000000 0.118716
H 0.000000 -0.756998 -0.474863


