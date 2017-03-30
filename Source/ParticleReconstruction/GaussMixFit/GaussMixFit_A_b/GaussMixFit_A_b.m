function [A b A_std Region]=GaussMixFit_A_b(x_Init, y_Init, A_Init, R_Init, b_Init, Range_A, Range_b, ...                                                    
                                            OtherParam, I, ObjLabel, Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior)

 
if isempty(OtherParam)
    MaxWindowRadius=9;
else
    MaxWindowRadius=OtherParam.MaxWindowRadius;
end

[Ly Lx]=size(I);

x=x_Init(:);
y=y_Init(:);
A=A_Init(:);
R=R_Init(:);
b=b_Init(:);

A_std=zeros(size(x));
%ParticleNum=length(x);

%% Select the region that contains the particles 

Region=SelectGaussMixRegion(I, ObjLabel, x, y, A, R, b, MaxWindowRadius);

RegionNum=length(Region);

ImageInfo.Lx=Lx;
ImageInfo.Ly=Ly;

for n=1:RegionNum
    
    %ParticleNum_n=length(Region(n).x);    
    ParticleIdx_n=Region(n).ParticleIdx(:);
        
    PixelIdxList_n=Region(n).PixelIdxList;
    
    RegionInfo.PixelIdxList=Region(n).PixelIdxList;
    RegionInfo.XPixelList=Region(n).XPixelList;
    RegionInfo.YPixelList=Region(n).YPixelList;

    %% Initialize the variables        
    x_n=x(ParticleIdx_n);
    y_n=y(ParticleIdx_n);
    A_n=A(ParticleIdx_n);
    R_n=R(ParticleIdx_n);
    b_n=b(ParticleIdx_n);         
                        
    Range_A_n=Range_A(:, ParticleIdx_n);
    Range_b_n=Range_b(:, ParticleIdx_n);
    
    I_Data=I(PixelIdxList_n);
    %% Estimation   
    if isempty(Ib_Prior) && isempty(A_Prior)
        [A_n b_n beta_n] = GaussMixFit_Region_A_b_NoPrior(x_n, y_n, A_n, R_n, b_n, Range_A_n, Range_b_n, ...
                                                          ImageInfo, RegionInfo, I_Data);    
    elseif ~isempty(Ib_Prior) && isempty(A_Prior)
        Ib_Prior_n=Ib_Prior(PixelIdxList);    
        Ib_SigmaPrior_n=Ib_SigmaPrior(PixelIdxList);
        [A_n b_n beta_n] = GaussMixFit_Region_A_b_with_b_Prior(x_n, y_n, A_n, R_n, b_n, Range_A_n, Range_b_n, ...
                                                               ImageInfo, RegionInfo, I_Data, Ib_Prior_n, Ib_SigmaPrior_n);
    elseif ~isempty(Ib_Prior) && ~isempty(A_Prior)                  
        
        Ib_Prior_n=Ib_Prior(PixelIdxList);    
        Ib_SigmaPrior_n=Ib_SigmaPrior(PixelIdxList);
        
        A_Prior_n=A_Prior(ParticleIdx_n);
        A_SigmaPrior_n=A_SigmaPrior(ParticleIdx_n);
            
        [A_n b_n beta_n] = GaussMixFit_Region_A_b_with_A_b_Prior(x_n, y_n, A_n, R_n, b_n, Range_A_n, Range_b_n, ...
                                                                 ImageInfo, RegionInfo, I_Data, Ib_Prior_n, Ib_SigmaPrior_n, A_Prior_n, A_SigmaPrior_n);
    end    
     
    [A_std_n MSE_n]=GaussMixFit_A_b_STATISTICS(x_n, y_n, A_n, R_n, beta_n, I_Data, RegionInfo.XPixelList, RegionInfo.YPixelList);
    
    %% Output x_n y_n a_n            
    x(ParticleIdx_n)=x_n;
    y(ParticleIdx_n)=y_n;
    A(ParticleIdx_n)=A_n;
    R(ParticleIdx_n)=R_n;
    b(ParticleIdx_n)=b_n;       
    
    A_std(ParticleIdx_n)=A_std_n;
    %% Update Rregion Information        
    Region(n).x=x_n;
    Region(n).y=y_n;        
    Region(n).A=A_n;
    Region(n).R=R_n;
    Region(n).b=b_n;
    
    Region(n).A_std=A_std_n;
    Region(n).MSE=MSE_n;
    
end