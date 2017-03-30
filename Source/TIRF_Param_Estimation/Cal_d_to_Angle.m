function Angle=Cal_d_to_Angle(d, nt, ni, Lambda)

Angle=asin(sqrt(((Lambda./(4*pi*d)).^2+nt^2)/ni^2));