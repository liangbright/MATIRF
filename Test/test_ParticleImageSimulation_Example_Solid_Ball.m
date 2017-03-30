Depth=1000;
AngleNum=length(Depth);
ni=1.522;
nt=1.333;
Lambda=488;
I0=Cal_d_to_I0(Depth, nt, ni, Lambda);
I0=I0/I0(1);
NA=1.49;

SysGain=50;

PixelSize=160;
TypicalParticleRadius=50;

MATIRF_Param.LateralResolution=0.61*Lambda/NA;
MATIRF_Param.PSFSigmaXY=170;
MATIRF_Param.PSFSigmaZ=inf;
MATIRF_Param.AngleNum=AngleNum;
MATIRF_Param.Depth=Depth;
MATIRF_Param.ni=ni;
MATIRF_Param.nt=nt;
MATIRF_Param.I0=I0;
MATIRF_Param.Lambda=Lambda;
MATIRF_Param.NA=NA;
MATIRF_Param.SysGain=SysGain;
MATIRF_Param.PixelSize=PixelSize;
%-----------------------------------------------------------------------------------------------------------
ParticleNum=4;
Lx=20; 
Ly=20;

x=[5;  5; 15; 15];
y=[5; 15; 15; 5];

R=[50; 50; 50; 50];

z=[50; 50; 50; 50];

C=[500; 666; 832; 1000];

CellThickness=1000;

ParticleFeature=zeros(5, ParticleNum);
ParticleFeature(1,:)=x*PixelSize;
ParticleFeature(2,:)=y*PixelSize;
ParticleFeature(3,:)=z;
ParticleFeature(4,:)=R;
ParticleFeature(5,:)=C;
%% Generate MATIRF Images
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
%
ImageSize=[Ly, Lx];

GridNumPerPixel=11;

tic
[Iobj, GridSize]=MATIRF_ParticleSimulation_numerical(MATIRF_Param, ParticleFeature, ImageSize, GridNumPerPixel);
toc
%%
BgDensity=10;
[ParticleSNRList, SNRPerAngle]=MATIRF_Cal_ParticleSNR(C, z, R, BgDensity, CellThickness, MATIRF_Param);

BgList=MATIRF_BackgroundSimulation(MATIRF_Param, CellThickness, BgDensity);
%%
imtool close all

bg=100;
I=Iobj+bg;
I(end,end)=10;
imtool(I)
%
I=poissrnd(Iobj+bg);
imtool(I)
I=I+sqrt(I).*randn(size(I)); 

I=I*100/bg;

imtool(I)