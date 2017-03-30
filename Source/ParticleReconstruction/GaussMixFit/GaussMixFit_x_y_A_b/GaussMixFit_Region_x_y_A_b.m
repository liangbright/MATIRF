function [x y A b Param] = GaussMixFit_Region_x_y_A_b(x_Init, y_Init, A_Init, R_Init, b_Init, Radius_xy, Range_A, Range_b, ...
                                                      ImageInfo, RegionInfo, I_Data, Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior)


Options = optimset('Jacobian','on','Display','off');

Lx=ImageInfo.Lx;
Ly=ImageInfo.Ly;

subX=RegionInfo.XPixelList;
subY=RegionInfo.YPixelList;  
    
% Initialize the variables
ParticleNum=length(x_Init);

Param=zeros(ParticleNum*5,1);

x=x_Init(:);
y=y_Init(:);
A=A_Init(:);
R=R_Init(:);
beta=b_Init(:);
          
x_min=max(x-Radius_xy(:), 0); x_max=min(x+Radius_xy(:), Lx);    
y_min=max(y-Radius_xy(:), 0); y_max=min(y+Radius_xy(:), Ly);  

A_min=Range_A(1,:)';
A_max=Range_A(2,:)';

beta_min=Range_b(1,:)';
beta_max=Range_b(2,:)';

Param_Init =[x;     y;     A;     beta];
Param_Lowbd=[x_min; y_min; A_min; beta_min];
Param_Upbd =[x_max; y_max; A_max; beta_max];
   
%% Estimation   
if nargin == 11
    
    Param=lsqnonlin(@GaussMixFit_x_y_A_b_FitErr_NoPrior, Param_Init, Param_Lowbd, Param_Upbd, Options, ...
                    ParticleNum, R, I_Data, subX, subY);    
elseif nargin == 13                      

    Param=lsqnonlin(@GaussMixFit_x_y_A_b_FitErr_with_b_Prior, Param_Init, Param_Lowbd, Param_Upbd, Options, ...
                    ParticleNum, R, I_Data, subX, subY, Ib_Prior, Ib_SigmaPrior);       
elseif nargin == 15

    Param=lsqnonlin(@GaussMixFit_x_y_A_b_FitErr_with_A_b_Prior, Param_Init, Param_Lowbd, Param_Upbd, Options, ...
                    ParticleNum, R, I_Data, subX, subY, Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior);   
end
    
%% Output x y A R b

sIdx=1; eIdx=ParticleNum;
x=Param(sIdx:eIdx);

sIdx=ParticleNum+1; eIdx=sIdx+ParticleNum-1;
y=Param(sIdx:eIdx);

sIdx=2*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
A=Param(sIdx:eIdx);

sIdx=3*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
beta=Param(sIdx:eIdx);
%%
b=Cal_b(x, y, A, R, beta, subX, subY);