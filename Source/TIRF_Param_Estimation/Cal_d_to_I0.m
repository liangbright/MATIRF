function I0=Cal_d_to_I0(d, nt, ni, Lambda)

Angle=Cal_d_to_Angle(d, nt, ni, Lambda);

I0=Cal_TIRF_I0(nt, ni, Angle);