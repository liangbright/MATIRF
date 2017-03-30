function [ImageStack, TotalIntensity, A]=MATIRF_ParticleSimulation(MATIRF_Param, ParticleFeature, ImageSize)

AngleNum=MATIRF_Param.AngleNum;
Depth=MATIRF_Param.Depth;
I0=MATIRF_Param.I0;

PixelSize=MATIRF_Param.PixelSize;

PSFSigmaXY=MATIRF_Param.PSFSigmaXY;

Ly=ImageSize(1);
Lx=ImageSize(2);

ImageStack=zeros(Ly, Lx, AngleNum);

x=ParticleFeature(1,:)/PixelSize;
y=ParticleFeature(2,:)/PixelSize;
z=ParticleFeature(3,:);
R=ParticleFeature(4,:)/PixelSize;
C=ParticleFeature(5,:);

ParticleNum=length(x);

A=zeros(ParticleNum, AngleNum);

TotalIntensity=zeros(ParticleNum, AngleNum);

tempR=sqrt(R.^2+PSFSigmaXY^2)/PixelSize;
%---------------------------------------------------------------------
for m=1:AngleNum                          
    
    d=Depth(m);               
            
    Alpha=3*((R/d).*cosh(R/d)-sinh(R/d)).*(d./R).^3;     

    Am=Alpha.*C.*I0(m).*exp(-z/d);       
               
    TotalIntensity(:,m)=Am*2*pi*(PSFSigmaXY/PixelSize)^2;
    
    I=GaussMix(x, y, Am, tempR, Lx, Ly);    
    
    ImageStack(:,:,m)=I;               
end