function [x y a r b] = GaussMixFit_Region_All(x_Init, y_Init, a_Init, r_Init, b_Init, Radius_xy, Range_r, ...
                                          ImageInfo, RegionInfo, I, Ib, Ib_Sigma)


Options = optimset('Jacobian','on','Display','off');

Lx=ImageInfo.Lx;
Ly=ImageInfo.Ly;

subX=RegionInfo.XPixelList;
subY=RegionInfo.YPixelList;  
PixelIdxList=RegionInfo.PixelIdxList;
    
% Initialize the variables
ParticleNum=length(x_Init);

Param=zeros(ParticleNum*5,1);

x=x_Init(:);
y=y_Init(:);
a=a_Init(:);
r=r_Init(:);
beta=b_Init(:);
          
x_min=max(x-Radius_xy, 0); x_max=min(x+Radius_xy, Lx);    
y_min=max(y-Radius_xy, 0); y_max=min(y+Radius_xy, Ly);  

a_min=zeros(ParticleNum, 1);
a_max=5*a+median(a);

r_min=Range_r(1)*ones(ParticleNum,1);
r_max=Range_r(2)*ones(ParticleNum,1);

beta_min=zeros(ParticleNum, 1);
beta_max=5*beta+median(beta);     
if nargin == 12   
    for k=1:ParticleNum
        yk=y(k); yk=round(yk); yk=min(max(yk,1), Ly);
        xk=x(k); xk=round(xk); xk=min(max(xk,1), Lx);
        bk=Ib(k);
        bk_std=Ib_Sigma(yk,xk);
        
        beta_min(k)=max(bk-5*bk_std, 0);
        beta_max(k)=min(bk+5*bk_std, beta_max(k));
    end
end

Param_Init =[x;     y;     a;     r;     beta];
Param_Lowbd=[x_min; y_min; a_min; r_min; beta_min];
Param_Upbd =[x_max; y_max; a_max; r_max; beta_max];
    
I_Data=I(PixelIdxList);    
%% Estimation   
if nargin == 10    
    
    Param=lsqnonlin(@GaussMix_FitErr_NoPrior, Param_Init, Param_Lowbd, Param_Upbd, Options, ...
                    ParticleNum, I_Data, subX, subY);    
elseif nargin == 12                      
    Ib_Prior=Ib(PixelIdxList);    
    Ib_SigmaPrior=Ib_Sigma(PixelIdxList);
    Param=lsqnonlin(@GaussMix_FitErr_with_b_Prior, Param_Init, Param_Lowbd, Param_Upbd, Options, ...
                    ParticleNum, I_Data, subX, subY, Ib_Prior, Ib_SigmaPrior);       
end
    
%% Output x_n y_n a_n    

sIdx=1; eIdx=ParticleNum;
x=Param(sIdx:eIdx);

sIdx=ParticleNum+1; eIdx=sIdx+ParticleNum-1;
y=Param(sIdx:eIdx);

sIdx=2*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
a=Param(sIdx:eIdx);

sIdx=3*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
r=Param(sIdx:eIdx);

sIdx=4*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
beta=Param(sIdx:eIdx);

b=Cal_b(x, y, a, r, beta, subX, subY);