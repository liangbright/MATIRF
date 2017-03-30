Depth=200;
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
MATIRF_Param.PSFSigmaXY=sqrt(PixelSize^2-TypicalParticleRadius^2);
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
Lx=100; 
Ly=100;

[x, y]=meshgrid(10:10:90, 10:10:90);
x=x(:);
y=y(:);
ParticleNum=length(x);
x=x+0.5*rand(ParticleNum, 1);
y=y+0.5*rand(ParticleNum, 1);
ParticleNum=length(x);
R=50*ones(ParticleNum, 1);

z=100*ones(ParticleNum, 1);
C=1000*ones(ParticleNum, 1);

CellThickness=2000;

ParticleFeature_outer=zeros(5, ParticleNum);
ParticleFeature_outer(1,:)=x*PixelSize;
ParticleFeature_outer(2,:)=y*PixelSize;
ParticleFeature_outer(3,:)=z;
ParticleFeature_outer(4,:)=R;
ParticleFeature_outer(5,:)=2*C;
%% Generate MATIRF Images
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
%%
ParticleFeature_inner=ParticleFeature_outer;
ParticleFeature_inner(4,:)=R-22;
ParticleFeature_inner(5,:)=C;
%%
ImageSize=[Ly, Lx];

for GridNumPerPixel=21

tic
[Iobj_outer, GridSize_outer]=MATIRF_ParticleSimulation_numerical(MATIRF_Param, ParticleFeature_outer, ImageSize, GridNumPerPixel);

[Iobj_inner, GridSize_inner]=MATIRF_ParticleSimulation_numerical(MATIRF_Param, ParticleFeature_inner, ImageSize, GridNumPerPixel);
toc
end
Iojb=Iobj_outer-Iobj_inner;
%%
imtool(Iojb)
%%
bg=500;
I=Iojb+bg;
I(end,end)=bg;
%imtool(I)
%
I=poissrnd(I);
%imtool(I)
I=I+sqrt(I).*randn(size(I)); 
%
imtool(I)