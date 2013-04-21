function [nonzeros res_forBit]=testomp;
!./ompbit
load Res_forBit0.mat
%load Res_forBit1.mat
%load Res_forBit2.mat
%load Res_forBit3.mat
%res_forBit = res_forBit0 + res_forBit1 + res_forBit2 + res_forBit3;
res_forBit = res_forBit0;
res_forBit = res_forBit' + res_forBit;
nonzeros = nnz(res_forBit)
