function [C z]=Estimate_InitParam(A, I0, Depth)

x=-1./Depth;

Y=log(A./I0);

X = [ones(size(x)) x];
Param = regress(Y, X);

C=exp(Param(1));
z=Param(2);