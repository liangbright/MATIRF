function ImageStack=MATIRF_NoiseSimulation(CleanImageStack, BgList, SysGain)

[Ly,Lx,AngleNum]=size(CleanImageStack);

ImageStack=zeros(Ly,Lx,AngleNum);

for m=1:AngleNum                                                                                 

    I=CleanImageStack(:,:,m);

    I=I+BgList(m);
    I=poissrnd(I);
    I=I+sqrt(I).*randn(Ly,Lx); 
    I=round(I*SysGain)/SysGain;
    
    ImageStack(:,:,m)=I;
end