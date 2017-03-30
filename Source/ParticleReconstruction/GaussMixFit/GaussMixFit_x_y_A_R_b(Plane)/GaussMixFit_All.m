function [x y a r b Region] = GaussMixFit_All(x_Init, y_Init, a_Init, r_Init, b_Init, SearchRadius_xy, Range_r, OtherParam, ...
                                          I, CellLabel, Ib, Ib_Sigma)


if isempty(OtherParam)
    MaxWindowRadius=9;
else
    MaxWindowRadius=OtherParam.MaxWindowRadius;
end

[Ly Lx]=size(I);

x=x_Init(:);
y=y_Init(:);
a=a_Init(:);
r=r_Init(:);
b=b_Init(:);
Radius_xy=SearchRadius_xy(:);

%ParticleNum=length(x);

%% Select the region that contains the particles 

Region=SelectGaussMixRegion(I, CellLabel, x, y, a, r, b, MaxWindowRadius);

RegionNum=length(Region);

ImageInfo.Lx=Lx;
ImageInfo.Ly=Ly;

for n=1:RegionNum
    
    %ParticleNum_n=length(Region(n).x);
    %PixelIdxList_n=Region(n).PixelIdxList;
    ParticleIdx_n=Region(n).ParticleIdx(:);
        
    RegionInfo.PixelIdxList=Region(n).PixelIdxList;
    RegionInfo.XPixelList=Region(n).XPixelList;
    RegionInfo.YPixelList=Region(n).YPixelList;

    %% Initialize the variables        
    x_n=x(ParticleIdx_n);
    y_n=y(ParticleIdx_n);
    a_n=a(ParticleIdx_n);
    r_n=r(ParticleIdx_n);
    b_n=b(ParticleIdx_n); 
    
    Radius_xy_n=Radius_xy(ParticleIdx_n);                      
    %% Estimation   
    if nargin == 10                      
        [x_n y_n a_n r_n b_n] = GaussMixFit_Region_All(x_n, y_n, a_n, r_n, b_n, Radius_xy_n, Range_r, ...
                                                   ImageInfo, RegionInfo, I);    
    elseif nargin == 12                        
        [x_n y_n a_n r_n b_n] = GaussMixFit_Region_All(x_n, y_n, a_n, r_n, b_n, Radius_xy_n, Range_r, ...
                                                   ImageInfo, RegionInfo, I, Ib, Ib_Sigma);    
    end
    
    %% Output x_n y_n a_n            
    x(ParticleIdx_n)=x_n;
    y(ParticleIdx_n)=y_n;
    a(ParticleIdx_n)=a_n;
    r(ParticleIdx_n)=r_n;
    b(ParticleIdx_n)=b_n;       
    
    %% Update Rregion Information        
    Region(n).x=x_n;
    Region(n).y=y_n;        
    Region(n).a=a_n;
    Region(n).r=r_n;
    Region(n).b=b_n;
    
end