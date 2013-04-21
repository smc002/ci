function testomp;
!./ompbit
load Res_forBit0.mat
load Res_forBit1.mat
load Res_forBit2.mat
load Res_forBit3.mat
res=res0+res1+res2+res3;
load Res_forBit.mat
%res_forBit = res_forBit0 + res_forBit1 + res_forBit2 + res_forBit3;
nnz(res_forBit-res)
