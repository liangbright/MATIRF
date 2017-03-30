function I0=Cal_TIRF_I0(nt, ni, Angle)
% p-polarization
nti=nt/ni;

temp1=4*cos(Angle).^2.*(2*sin(Angle).^2-nti^2);

temp2=nti^4*cos(Angle).^2+sin(Angle).^2-nti^2;

I0=temp1./temp2;