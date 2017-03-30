function Hfilter= Mask_TopHat(ParticleRadius)

R=ceil(4*ParticleRadius);

[X Y]=meshgrid(-R:R, -R:R);
dist=sqrt(X.^2+Y.^2);

f1=zeros(2*R+1, 2*R+1);
f2=zeros(2*R+1, 2*R+1);

f1(dist<=(2*ParticleRadius))=1;
f2(dist>(2*ParticleRadius))=-1;

f1=f1/sum(f1(:));
f2=f2/(sum(-f2(:)));

Hfilter= 2.4385*(f1+f2);