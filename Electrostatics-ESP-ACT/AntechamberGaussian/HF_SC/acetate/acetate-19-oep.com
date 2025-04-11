%mem=12000MB
%nprocshared=8
%chk=acetate19.chk
#P HF/aug-cc-pvtz Opt=(Redundant, calcall, verytight) symm=(loose,follow) Pop=(MK,Hirshfeld,ReadRadii) iop(6/33=2) iop(6/42=6) Polar Freq
maxdisk=128GB

acetate

-1 1
C -1.355140 -0.036667 -0.000539
C 0.221434 -0.000902 0.000025
O 0.790199 -1.115292 -0.012178
O 0.711240 1.151742 0.012731
H -1.740335 -1.061211 -0.011487
H -1.734589 0.488177 0.885015
H -1.734354 0.506852 -0.874862


