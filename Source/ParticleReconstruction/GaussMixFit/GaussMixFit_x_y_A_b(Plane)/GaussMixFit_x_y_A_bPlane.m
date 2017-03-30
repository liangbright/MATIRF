function [Param A x y b Ic RegionInfo Region F J] = GaussMix_Fit(PSFSigma, A_Init, x_Init, y_Init, b_Init, I, ...
                                                                 Ib, Ib_Sigma, A_Input, A_Sigma)

[Ly Lx]=size(I);

A=A_Init(:);
x=x_Init(:);
y=y_Init(:);

ParticleNum=length(x);

Options = optimset('Jacobian','on','Display','off');

%% Select the region that contains the particles    
[Region RegionInfo]=SelectGaussMixRegion(PSFSigma, A, x, y, Lx, Ly);
RegionNum=length(Region);

Param=zeros(ParticleNum*6,1);

b_Init_Median=median(b_Init);

A_Init_Median=median(A_Init);

for n=1:RegionNum
    
    ParticleNum_n=length(Region(n).x);
    PixelIdxList=Region(n).PixelIdxList;
    ParticleIdx=Region(n).ParticleIdx(:);
    
    %% Initialize the variables    
    A_n=Region(n).A;
    x_n=Region(n).x;
    y_n=Region(n).y;        
    b_n=b_Init(ParticleIdx); 
    
    I_Data=I(PixelIdxList);
    subX=Region(n).XPixelList;
    subY=Region(n).YPixelList;                
    
    P0_Init=b_n;
    P0_Min=zeros(size(P0_Init));
    P0_Max=2*b_n+b_Init_Median;
    
    P_Init=0.1*ones(ParticleNum_n,1);

    const=PSFSigma+0.5;
    x_n_min=max(x_n-const, 0); x_n_max=min(x_n+const, Lx);    
    y_n_min=max(y_n-const, 0); y_n_max=min(y_n+const, Ly);    
    
    Param_Init =[A_n(:);                  P0_Init;  0*P_Init; 0*P_Init; x_n;         y_n];
    Param_Lowbd=[0*A_n(:);                 P0_Min;   -P_Init;  -P_Init; x_n_min; y_n_min];
    Param_Upbd =[5*A_n(:)+A_Init_Median;   P0_Max;    P_Init;   P_Init; x_n_max; y_n_max];

    %% Estimation   
    if nargin == 6
        Param_n=lsqnonlin(@GaussMix_FitErr_NoPrior, Param_Init, Param_Lowbd, Param_Upbd, Options, ...
                          ParticleNum_n, PSFSigma, I_Data, subX, subY);
    elseif nargin == 8        
        Ib_Prior=Ib(PixelIdxList);
        Ib_SigmaPrior=Ib_Sigma(PixelIdxList);   
        
        Param_n=lsqnonlin(@GaussMix_FitErr_with_b_Prior, Param_Init, Param_Lowbd, Param_Upbd, Options, ...
                          ParticleNum_n, PSFSigma, I_Data, subX, subY, Ib_Prior, Ib_SigmaPrior);
    elseif nargin == 10
        Ib_Prior=Ib(PixelIdxList);
        Ib_SigmaPrior=Ib_Sigma(PixelIdxList);        
        A_Prior=A_Input(ParticleIdx);
        A_SigmaPrior=A_Sigma(ParticleIdx);
        
        Param_n=lsqnonlin(@GaussMix_FitErr_with_A_b_Prior, Param_Init, Param_Lowbd, Param_Upbd, Options, ...
                          ParticleNum_n, PSFSigma, I_Data, subX, subY, Ib_Prior, Ib_SigmaPrior, A_Prior, A_SigmaPrior);                
    end
    
    %% Output x_n y_n A_n    
    [A_n x_n y_n P0_n P1_n P2_n]=SeperateParam(Param_n, ParticleNum_n, 1);        
    
    Param=SetParam(Param, ParticleNum, 1, Region(n).ParticleIdx, A_n, x_n, y_n, P0_n, P1_n, P2_n);    
end

%% Output x y A
[A x y]=SeperateParam(Param, ParticleNum, 1); 

%%
subX=RegionInfo.XPixelList;
subY=RegionInfo.YPixelList;
subIdxList=RegionInfo.PixelIdxList;

I_Data=I(subIdxList);

[F J b_Model]=GaussMix_FitErr(Param, ParticleNum, PSFSigma, I_Data, subX, subY);

%%
b=zeros(ParticleNum, 1);
Ic=I;
Ic(subIdxList)=b_Model;
for k=1:ParticleNum
    yk=y(k); yk=round(yk); yk=min(max(yk,1), Ly);
    xk=x(k); xk=round(xk); xk=min(max(xk,1), Lx);
    b(k)=Ic(yk, xk);
end
%%