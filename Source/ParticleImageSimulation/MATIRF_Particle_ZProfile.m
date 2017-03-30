function Value=MATIRF_Particle_ZProfile(z, x, y, MATIRF_Param, AngleIndex, ParticleFeature)

Depth=MATIRF_Param.Depth;

I0=MATIRF_Param.I0;

x0=ParticleFeature(1);
y0=ParticleFeature(2);
z0=ParticleFeature(3);
R=ParticleFeature(4);
C=ParticleFeature(5);

Volume=(4*pi/3)*R.^3;

dist_sq=(x-x0).^2+(y-y0).^2+(z-z0).^2;

Value=(C/Volume)*I0(AngleIndex)*exp(-z/Depth(AngleIndex));

Value(dist_sq>R^2)=0;