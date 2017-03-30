function Hfilter= Mask_Gauss(ParticleRadius)

R=ceil(4*ParticleRadius);

%Sigma=Sigma*0.6745;

[X Y]=meshgrid(-R:R, -R:R);
temp=X.^2+Y.^2;

Hfilter=exp(-temp/(2*ParticleRadius^2));

Hfilter=Hfilter/sum(Hfilter(:));