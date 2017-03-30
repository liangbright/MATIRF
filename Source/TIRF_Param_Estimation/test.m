ni=1.522;
nt=1.333;
Lambda=488;

d=[64 100 200 300 400];

Angle=Cal_d_to_Angle(d, nt, ni, Lambda);

Angle=[Angle 1];

I0=Cal_TIRF_I0(nt, ni, Angle);
figure; plot(I0)

%I0=I0/I0(1)