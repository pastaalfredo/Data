%mem=12000MB
%nprocshared=8
%chk=ethylammonium19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB

ethylammonium

1 1
N -1.204745 -0.272143 0.011120
C 0.047392 0.609966 0.009339
H -2.064013 0.285530 0.031323
H -1.216615 -0.892820 0.826661
H -1.237932 -0.865880 -0.823661
H -0.021991 1.223989 0.907303
H -0.044888 1.252900 -0.866044
C 1.303284 -0.240324 -0.020751
H 1.383371 -0.883693 0.859113
H 2.170935 0.421837 -0.021220
H 1.360291 -0.854714 -0.922844


