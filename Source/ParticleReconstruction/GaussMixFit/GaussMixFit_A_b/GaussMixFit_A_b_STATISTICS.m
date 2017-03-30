function [A_std MSE]=GaussMixFit_A_b_STATISTICS(x, y, A, R, beta, I_Data, XPixelList, YPixelList)

ParticleNum=length(x);
PixelNum=length(XPixelList);

%--------------------------------------------------------------------------------------
Coef=zeros(PixelNum, 2*ParticleNum);
Iobj=zeros(PixelNum,1);
Ib=zeros(PixelNum,1);
b_Weight=zeros(PixelNum,ParticleNum); 
    
for k=1:ParticleNum    
    tempEXP=exp(-((x(k)-XPixelList).^2+(y(k)-YPixelList).^2)/(2*R(k)^2));  
    Coef(:,k)=tempEXP;
    Iobj=Iobj+A(k)*tempEXP;
end
               
for k=1:ParticleNum                
    b_Weight(:,k)=A(k)*Coef(:,k)+eps;  
    Ib=Ib+beta(k).*b_Weight(:,k);
end    
tempSum=sum(b_Weight,2);
Ib=Ib./tempSum;
    
for n=1:PixelNum    
    b_Weight(n,:)=b_Weight(n,:)/tempSum(n);
end
    
Coef(:,ParticleNum+1:2*ParticleNum)=b_Weight;
    
tempVar=2*(Iobj+Ib);
WI=1./sqrt(tempVar);    
    
for k=1:2*ParticleNum
    Coef(:,k)=Coef(1:PixelNum,k).*WI;          
end    
    
Data=I_Data.*WI;

%--------------------------------------------------------------------------------------------
A_std=zeros(ParticleNum, 1);
%b_std=zeros(ParticleNum, 1);

IdxListP=find(A>0.001);
IdxListN=find(A<=0.001);

if ~isempty(IdxListP)
    IdxList=[IdxListP; ParticleNum+IdxListP];

    F=Coef(:,IdxList)*[A(IdxListP); beta(IdxListP)]-Data;

    NumDegFree=PixelNum-2*length(IdxList);
    Param_Std=sqrt(diag(inv(Coef(:,IdxList)'*Coef(:,IdxList))*sum(F.^2)/NumDegFree));

    NumP=length(IdxListP);
    A_std(IdxListP)=Param_Std(1:NumP);   
    A_std(IdxListN)=max(A_std(IdxListP));
else
    A_std(:)=NaN;
end
%b_std(A_IdxListP)=Param_Std(NumP+1:2*NumP);
%b_std(A_IdxListN)=max(Param_Std(NumP+1:2*NumP));
%--------------------------------------------------------------------------------------------
MSE=mean((Iobj+Ib-I_Data).^2);
