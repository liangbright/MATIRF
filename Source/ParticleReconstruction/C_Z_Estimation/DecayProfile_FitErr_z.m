function [F J]=DecayProfile_FitErr_z(Param, C, A, I0, Depth, W)

z=Param;
tempEXP=exp(-z./Depth);
F=(I0.*C.*tempEXP-A).*W;
J=I0.*C.*tempEXP.*(-1./Depth).*W;