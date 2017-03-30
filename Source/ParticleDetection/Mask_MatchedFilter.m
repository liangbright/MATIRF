function Hfilter= Mask_MatchedFilter(ParticleRadius)

R=ceil(4*ParticleRadius);

[X Y]=meshgrid(-R:R, -R:R);
distMap=X.^2+Y.^2;

f1=(1-distMap/((2*ParticleRadius)^2)).*exp(-distMap/(2*ParticleRadius^2));

fn_1=f1;
fn_1(f1>0)=0;
fn_1=fn_1/sum(fn_1(:));

f2=exp(-distMap/(2*ParticleRadius^2));

fp_2=f2;
fp_2(f1<0)=0;
fp_2=fp_2/sum(fp_2(:));

Hfilter=1.8802*(fp_2-fn_1);