function [A b beta]=GaussMixFit_Region_A_b_with_b_Prior(x_Init, y_Init, A_Init, R_Init, b_Init, Range_A, Range_b, ...                                                    
                                                        ImageInfo, RegionInfo, I_Data, Ib_Prior, Ib_SigmaPrior)

ParticleNum=length(x_Init);

x=x_Init(:);
y=y_Init(:);
A=A_Init(:);
R=R_Init(:);
beta=b_Init(:);

Amin=Range_A(1,:)';
Amax=Range_A(2,:)';

beta_min=Range_b(1,:)';
beta_max=Range_b(2,:)';

subIdxList=RegionInfo.PixelIdxList;
subX=RegionInfo.XPixelList;subY=RegionInfo.YPixelList;
    
PixelNum=length(subIdxList);   
    
MaxIteration=2;
for Iteration=1:MaxIteration
        
    if Iteration < MaxIteration
        Options = optimset('MaxIter', 100, 'Display','off');
    else
        Options = optimset('MaxIter', 200, 'Display','off');
    end
    
    %% Initialize the variables               
    Coef=zeros(PixelNum+ParticleNum, 2*ParticleNum);
    Iobj=zeros(PixelNum,1);
    Ib=zeros(PixelNum,1);
    b_Weight=zeros(PixelNum,ParticleNum); 
    
    for k=1:ParticleNum    
        tempEXP=exp(-((x(k)-subX).^2+(y(k)-subY).^2)/(2*R(k)^2));
        Coef(1:PixelNum, k)=tempEXP;
        Iobj=Iobj+A(k)*tempEXP;
    end
              
    for k=1:ParticleNum                
        b_Weight(:,k)=A(k)*Coef(1:PixelNum,k)+eps;  
        Ib=Ib+beta(k).*b_Weight(:,k);
    end        
    tempSum=sum(b_Weight,2);
    Ib=Ib./tempSum;
    
    for n=1:PixelNum    
        b_Weight(n,:)=b_Weight(n,:)/tempSum(n);
    end    
    Coef(1:PixelNum, ParticleNum+1:2*ParticleNum)=b_Weight;
    
    tempVar=2*(Iobj+Ib);
    WI=1./sqrt(tempVar+1);    
    
    Wb_Prior=1./(Ib_SigmaPrior+1);   
    
    %--------------------------------------------------------------------------------------------
    for k=1:2*ParticleNum
        Coef(1:PixelNum,k)=Coef(1:PixelNum,k).*WI;          
    end    
    Coef(PixelNum+1:PixelNum+ParticleNum, ParticleNum+1:2*ParticleNum)=diag(Wb_Prior);    
    
    Data=[I_Data.*WI; Ib_Prior.*Wb_Prior];
     
    %% ---------------------------------------------------------------------------------------    
    Param_Init=[A; beta];
    Lowbd=[Amin; beta_min];
    Upbd=[Amax;  beta_max];
    Param= lsqlin(Coef, Data, [], [], [], [], Lowbd, Upbd, Param_Init, Options);

    %%
    A=Param(1:ParticleNum);
    Amax=10*A+eps;
    Amin=zeros(ParticleNum, 1);

    beta=Param(ParticleNum+1:2*ParticleNum);    
    beta_max=10*beta+eps;
    beta_min=zeros(ParticleNum, 1);

end

b=Cal_b(x, y, A, R, beta, subX, subY);