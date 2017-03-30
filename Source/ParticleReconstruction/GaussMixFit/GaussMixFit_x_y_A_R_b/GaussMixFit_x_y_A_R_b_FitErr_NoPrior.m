function [F J Iobj Ib]=GaussMixFit_x_y_A_R_b_FitErr_NoPrior(Param, ParticleNum, I_Data, X, Y)
% Param = [x(1:ParticleNum) y(1:ParticleNum) a(1:ParticleNum) r(1:ParticleNum) beta(1:ParticleNum)]
% each particle has (x, y, a, r, b)

PixelNum=length(X);

sIdx=1; eIdx=ParticleNum;
x=Param(sIdx:eIdx);

sIdx=ParticleNum+1;   eIdx=sIdx+ParticleNum-1;
y=Param(sIdx:eIdx);

sIdx=2*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
A=Param(sIdx:eIdx);

sIdx=3*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
R=Param(sIdx:eIdx);

sIdx=4*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
beta=Param(sIdx:eIdx);

%---------------------------------------------------------------------------------------------
F=zeros(PixelNum,1);
J=zeros(PixelNum, 5*ParticleNum);
%%

Iobj=zeros(PixelNum,1);
Ib=zeros(PixelNum,1);
b_Weight=zeros(PixelNum,ParticleNum);

dFdx=zeros(PixelNum,ParticleNum);    
dFdy=zeros(PixelNum,ParticleNum);    
dFdA=zeros(PixelNum,ParticleNum);    
dFdR=zeros(PixelNum,ParticleNum);    
dFdb=zeros(PixelNum,ParticleNum);

for k=1:ParticleNum
    
    tempDist=(x(k)-X).^2+(y(k)-Y).^2;
    
    tempEXP=exp(-tempDist/(2*R(k)^2));      
        
    Iobj=Iobj+A(k)*tempEXP;
        
    b_Weight(:,k)=A(k)*tempEXP+eps;
    Ib=Ib+beta(k).*b_Weight(:,k);      
               
    dFdx(:,k)=A(k)*tempEXP.*((X-x(k))/(R(k)^2));
    
    dFdy(:,k)=A(k)*tempEXP.*((Y-y(k))/(R(k)^2));        
    
    dFdA(:,k)=tempEXP;
    
    dFdR(:,k)=A(k)*tempEXP.*tempDist/(R(k)^3);
        
end     

tempSum=sum(b_Weight,2);
for n=1:PixelNum    
    b_Weight(n,:)=b_Weight(n,:)/tempSum(n);
end
Ib=Ib./tempSum;

tempVar=2*(Iobj+Ib);
WI=1./sqrt(tempVar+1);
%WI=1;

for k=1:ParticleNum    
    
    dFdx(:,k)=dFdx(:,k).*WI;           
    dFdy(:,k)=dFdy(:,k).*WI;   
    dFdA(:,k)=dFdA(:,k).*WI;                
    dFdR(:,k)=dFdR(:,k).*WI;
                               
    dFdb(:,k)=b_Weight(:,k).*WI;                 
end

%--------------------------------------------------------------------------------
F=(Iobj+Ib-I_Data).*WI;
%--------------------------------------------------------------------------------
    
%--------------------------------------------------------------------------------
J=[dFdx dFdy dFdA dFdR dFdb]; 
%--------------------------------------------------------------------------------