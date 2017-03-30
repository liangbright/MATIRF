function [F J]=DecayProfile_FitErr(Param, A, I0, Depth, W)

C=Param(1);
z=Param(2);

tempEXP=exp(-z./Depth);
F=(I0.*C.*tempEXP-A).*W;

J=zeros(length(A),2);
J(:,1)=I0.*tempEXP.*W;
J(:,2)=I0.*C.*tempEXP.*(-1./Depth).*W;