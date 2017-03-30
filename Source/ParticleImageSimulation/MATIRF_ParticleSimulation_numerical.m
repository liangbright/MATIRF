function [ImageStack, TotalIntensity, GridSize, GridNumPerPixel]=MATIRF_ParticleSimulation_numerical(MATIRF_Param, ...
                                                                                     ParticleFeature, ImageSize, ...
                                                                                     GridNumPerPixel_Init)

PSFSigmaXY=MATIRF_Param.PSFSigmaXY;

AngleNum=MATIRF_Param.AngleNum;

PixelSize=MATIRF_Param.PixelSize;

Ly=ImageSize(1);
Lx=ImageSize(2);

ImageStack=zeros(Ly, Lx, AngleNum);
%------------------------------------------------------------------------------------------
x=ParticleFeature(1,:);
y=ParticleFeature(2,:);
z=ParticleFeature(3,:);
R=ParticleFeature(4,:);
C=ParticleFeature(5,:);

ParticleNum=length(x);
%------------------------------------------------------------------------------------------
% Scale : Control Simulation Precision
GridNumPerPixel=2*round(GridNumPerPixel_Init/2)-1;
GridSize=PixelSize/GridNumPerPixel;
%----------------------------------------------------------------------------------------------------
temp_PSFSigmaXY=PSFSigmaXY/GridSize;

W=ceil(4*temp_PSFSigmaXY);

[tempY, tempX]=ndgrid(-W:W, -W:W);

temp_xy=(tempX/temp_PSFSigmaXY).^2+(tempY/temp_PSFSigmaXY).^2;
PSFMask_xy=exp(-0.5*temp_xy);
%----------------------------------------------------------------------------------------------------
SumFunction=@(block_struct) sum(block_struct.data(:));
%------------------------------------------------------------------------------------------
TotalIntensity=zeros(ParticleNum,AngleNum);

for k=1:ParticleNum
    
    %delta=3*PSFSigmaXY+R(k);
    delta=3*PSFSigmaXY+2*R(k);
    
    temp_x1=round((x(k)-delta)/PixelSize); temp_x1=max(temp_x1, 1);
    temp_x2=round((x(k)+delta)/PixelSize); temp_x2=min(temp_x2, Lx);
    
    temp_y1=round((y(k)-delta)/PixelSize); temp_y1=max(temp_y1, 1);
    temp_y2=round((y(k)+delta)/PixelSize); temp_y2=min(temp_y2, Ly);    
    
    x1=temp_x1*PixelSize; y1=temp_y1*PixelSize;
    x2=temp_x2*PixelSize; y2=temp_y2*PixelSize;
    
    xList=(x1-GridSize*(GridNumPerPixel-1)/2):GridSize:(x2+GridSize*(GridNumPerPixel-1)/2);
    yList=(y1-GridSize*(GridNumPerPixel-1)/2):GridSize:(y2+GridSize*(GridNumPerPixel-1)/2);
              
    IdxList_x=find(xList>x(k)-R(k) & xList<x(k)+R(k));
    IdxList_y=find(yList>y(k)-R(k) & yList<y(k)+R(k));
    
    z_min=max(z(k)-R(k), 0);
    z_max=z(k)+R(k);
    
    tempParticleFeature=ParticleFeature(:,k);
    tempParticleFeature(5)=max(tempParticleFeature(5), 1);
    
    for m=1:AngleNum
        
        Im=MATIRF_ZProjection_2DImage(length(xList), length(yList), IdxList_x, IdxList_y, xList(IdxList_x), yList(IdxList_y),...
                                      z_min, z_max, ...
                                      tempParticleFeature, MATIRF_Param, m);
        
        tempTotalIntensity=MATIRF_Particle_3DIntegration(ParticleFeature(:,k), MATIRF_Param, m); 

       % imtool(Im)
        
        Im=imfilter(Im, PSFMask_xy);     

       % imtool(Im)

        if GridNumPerPixel > 1
            I=blockproc(Im, [GridNumPerPixel, GridNumPerPixel], SumFunction);
            I=I/(GridNumPerPixel)^2;
        else
            I=Im;
        end
   
        tempTotalIntensity=tempTotalIntensity*2*pi*(PSFSigmaXY/PixelSize)^2;        
        I=I*(tempTotalIntensity/sum(I(:)));
        
       % imtool(I)
        
        TotalIntensity(k,m)=tempTotalIntensity;
        
        ImageStack(temp_y1:temp_y2, temp_x1:temp_x2, m)=ImageStack(temp_y1:temp_y2, temp_x1:temp_x2, m)+I;
    end
end                                                   