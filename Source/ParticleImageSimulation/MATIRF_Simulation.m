function [ImageStack, A]=MATIRF_Simulation(MATIRF_Param, ParticleFeature, ImageSize, CellThickness, BgDensity)

AngleNum=MATIRF_Param.AngleNum;
Depth=MATIRF_Param.Depth;
I0=MATIRF_Param.I0;

PixelSize=MATIRF_Param.PixelSize;

PSFSigmaXY=MATIRF_Param.PSFSigmaXY;

SysGain=MATIRF_Param.SysGain;


Ly=ImageSize(1);
Lx=ImageSize(2);

ImageStack=zeros(Ly, Lx, AngleNum);

bg=zeros(AngleNum, 1);

x=ParticleFeature(1,:)/PixelSize;
y=ParticleFeature(2,:)/PixelSize;
z=ParticleFeature(3,:);
R=ParticleFeature(4,:)/PixelSize;
C=ParticleFeature(5,:);

ParticleNum=length(x);

A=zeros(ParticleNum, AngleNum);

tempR=sqrt(R.^2+PSFSigmaXY^2)/PixelSize;
%---------------------------------------------------------------------
for m=1:AngleNum                          
    
    d=Depth(m);               
            
    %Alpha=3*((R/d).*cosh(R/d)-sinh(R/d)).*(d./R).^3;     
    Alpha=1;
    
    Am=I0(m)*C.*Alpha.*exp(-z/d);       
    
    A(:,m)=Am(:);
        
    Iobj=GaussMix(x, y, Am, tempR, Lx, Ly);
        
    bg(m)=I0(m)*min(Depth(m), CellThickness)*BgDensity;
     
    %bg(m)=bg(m)*2*pi*(PSFSigmaXY/PixelSize)^2;
    
    Im=Iobj+bg(m);
    Im=poissrnd(Im);
    Im=Im+sqrt(Im).*randn(Ly,Lx); 
    Im=round(Im*SysGain)/SysGain;
    
    ImageStack(:,:,m)=Im;               
end