function d=Cal_TIRF_d(nt, ni, Angle, Lambda)

d=(Lambda/(4*pi))./sqrt(ni^2*sin(Angle).^2-nt^2);