function BgList=MATIRF_BackgroundSimulation(MATIRF_Param, CellThickness, BgDensity)

AngleNum=MATIRF_Param.AngleNum;
Depth=MATIRF_Param.Depth;
I0=MATIRF_Param.I0;
PSFSigmaXY=MATIRF_Param.PSFSigmaXY;
PixelSize=MATIRF_Param.PSFSigmaXY;

BgList=zeros(1,AngleNum);

for m=1:AngleNum
        
    BgList(m)=BgDensity*I0(m)*min(Depth(m), CellThickness);
    
    %{
    % integration from z=0 to z=CellThickness
    BgList(m)=I0(m)*min(Depth(m), CellThickness);
    % integration at each grid
    BgList(m)=BgList(m)*(BgDensity*PixelSize^2);    
    % PSF effect:
    BgList(m)=BgList(m)*2*pi*(PSFSigmaXY/PixelSize)^2;             
    %}
end