function [NX, P]=epdfhist(X, nbins, Xmin, Xmax)

if nargin == 2
    Xmax=max(X);
    Xmin=min(X);
end

step = (Xmax-Xmin)/(nbins-1);
binc = Xmin:step:Xmax;     
[N, NX] = hist(X, binc);
P=N/sum(N)/step;