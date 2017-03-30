function [A b A_std b_std]=GaussMix_Fit_A_b_NoPrior_TEST(PSFSigma, x, y, I)

[Ly Lx]=size(I);

ParticleNum=length(x);
A=zeros(ParticleNum, 1);
b=zeros(ParticleNum, 1);
for k=1:ParticleNum        
    yk=round(y(k));
    xk=round(x(k));    
    A(k)=I(yk, xk);
    b(k)=0.1*A(k);
end

MaxIteration=10;
for Iteration=1:MaxIteration
        
    if Iteration < MaxIteration
        Options = optimset('MaxIter', 100, 'Display','off');
    else
        Options = optimset('MaxIter', 200, 'Display','off');
    end
        
    A_Init=A;    
    b_Init=b;
    
    %% Select the region that contains the particles    
    [subX subY subIdxList]=SelectGaussMixRegion(PSFSigma, x, y, A_Init, Lx, Ly);
    I_Data=I(subIdxList);    
    PixelNum=length(subIdxList);
        
    %% Initialize the variables               
    Coef=zeros(PixelNum, 2*ParticleNum);
    Iobj=zeros(PixelNum,1);
    Ib=zeros(PixelNum,1);
    b_Weight=zeros(PixelNum,ParticleNum);  
    
    for k=1:ParticleNum    
        tempEXP=exp(-((x(k)-subX).^2+(y(k)-subY).^2)/(2*PSFSigma^2));
        Coef(:, k)=tempEXP;
        
        Iobj=Iobj+A(k)*tempEXP;
    end
                                  
    for k=1:ParticleNum                
        b_Weight(:,k)=A(k)*Coef(:,k)+eps;  
        Ib=Ib+b(k).*b_Weight(:,k);  
    end        
    tempSum=sum(b_Weight,2);
    Ib=Ib./tempSum;
    
    for n=1:PixelNum    
        b_Weight(n,:)=b_Weight(n,:)/tempSum(n);
    end        
    Coef(:, ParticleNum+1:2*ParticleNum)=b_Weight;
    
    tempVar=2*(Iobj+Ib);
    WI=1./sqrt(tempVar);               
    %WI=1-exp(-sqrt(Iobj));
    %WI=1-exp(-0.5*Iobj.^2./(2*(Ib)));
    %WI=sqrt(WI);
    %--------------------------------------------------------------------------------------------
    for k=1:2*ParticleNum
        Coef(:,k)=Coef(1:PixelNum,k).*WI;          
    end    
    
    Data=I_Data.*WI;
    
    %% ---------------------------------------------------------------------------------------           
    Param_Init=[A_Init; b_Init];
    Lowbd=0*Param_Init;
    Upbd=10*Param_Init;    
    Param= lsqlin(Coef, Data, [], [], [], [], Lowbd, Upbd, Param_Init, Options);
    A=Param(1:ParticleNum);
    b=Param(ParticleNum+1:2*ParticleNum);
end

A_std=zeros(ParticleNum, 1);
b_std=zeros(ParticleNum, 1);

A_IdxListP=find(A>0.001);
A_IdxListN=find(A<=0.001);

IdxList=[A_IdxListP; ParticleNum+A_IdxListP];

Residual=Coef(:,IdxList)*Param(IdxList)-Data;
NumDegFree=PixelNum-2*length(IdxList);
Param_Std=sqrt(diag(inv(Coef(:,IdxList)'*Coef(:,IdxList))*sum(Residual.^2)/NumDegFree));

NumP=length(A_IdxListP);
A_std(A_IdxListP)=Param_Std(1:NumP);
b_std(A_IdxListP)=Param_Std(NumP+1:2*NumP);

A_std(A_IdxListN)=mean(Param_Std(1:NumP));
b_std(A_IdxListN)=mean(Param_Std(NumP+1:2*NumP));