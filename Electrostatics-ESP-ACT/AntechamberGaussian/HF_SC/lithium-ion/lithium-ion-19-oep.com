%mem=12000MB
%nprocshared=8
%chk=lithium-ion19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB

lithium-ion

1 1
Li 0.000000 0.000000 0.000000


