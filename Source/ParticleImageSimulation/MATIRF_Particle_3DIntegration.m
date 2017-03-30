function Value=MATIRF_Particle_3DIntegration(ParticleFeature, MATIRF_Param, AngleIndex) 

d=MATIRF_Param.Depth(AngleIndex);
I0=MATIRF_Param.I0(AngleIndex);

z=ParticleFeature(3,:);
R=ParticleFeature(4,:);
C=ParticleFeature(5,:);

Value=3.*((R/d).*cosh(R/d)-sinh(R/d)).*(d./R).^3;
Value=Value.*C.*I0.*exp(-z/d);