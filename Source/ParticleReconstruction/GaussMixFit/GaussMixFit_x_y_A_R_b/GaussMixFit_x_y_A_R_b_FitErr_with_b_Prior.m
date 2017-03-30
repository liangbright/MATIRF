function [F J Iobj Ib]=GaussMixFit_x_y_A_R_b_FitErr_with_b_Prior(Param, ParticleNum, I_Data, X, Y, Ib_Prior, Ib_SigmaPrior)
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
F=zeros(2*PixelNum, 1);
J=zeros(2*PixelNum, 5*ParticleNum);
%%

Iobj=zeros(PixelNum,1);
Ib=zeros(PixelNum,1);
b_Weight=zeros(PixelNum,ParticleNum);

dFdx=zeros(2*PixelNum,ParticleNum);    
dFdy=zeros(2*PixelNum,ParticleNum);    
dFdA=zeros(2*PixelNum,ParticleNum);    
dFdR=zeros(2*PixelNum,ParticleNum);    
dFdb=zeros(2*PixelNum,ParticleNum);

for k=1:ParticleNum
    
    tempDist=(x(k)-X).^2+(y(k)-Y).^2;
    
    tempEXP=exp(-tempDist/(2*R(k)^2));       
        
    Iobj=Iobj+A(k)*tempEXP;
        
    b_Weight(:,k)=A(k)*tempEXP+eps;
    Ib=Ib+beta(k).*b_Weight(:,k);      
    
    dFdx(1:PixelNum,k)=A(k)*tempEXP.*((X-x(k))/(R(k)^2));
    
    dFdy(1:PixelNum,k)=A(k)*tempEXP.*((Y-y(k))/(R(k)^2)); 
        
    dFdA(1:PixelNum,k)=tempEXP;
    
    dFdR(1:PixelNum,k)=A(k)*tempEXP.*tempDist/(R(k)^3);        
    
end     

tempSum=sum(b_Weight,2);
for n=1:PixelNum    
    b_Weight(n,:)=b_Weight(n,:)/tempSum(n);
end
Ib=Ib./tempSum;

tempVar=2*(Iobj+Ib);
WI=1./sqrt(tempVar+1); 

WIb=1./(Ib_SigmaPrior+1);
    
for k=1:ParticleNum    
    
    dFdx(1:PixelNum,k)=dFdx(1:PixelNum,k).*WI;           
    dFdy(1:PixelNum,k)=dFdy(1:PixelNum,k).*WI;   
    dFdA(1:PixelNum,k)=dFdA(1:PixelNum,k).*WI;            
    dFdR(1:PixelNum,k)=dFdR(1:PixelNum,k).*WI;
                                
    dFdb(1:PixelNum,k)=b_Weight(:,k).*WI;
    
    % Cost=sum{(Ib-b_Data).^2}    
    sIdx=PixelNum+1; eIdx=2*PixelNum;
    dFdb(sIdx:eIdx,k)=b_Weight(:,k).*WIb;
end
    
%--------------------------------------------------------------------------------
F(1:PixelNum)=(Iobj+Ib-I_Data).*WI;
F(PixelNum+1:2*PixelNum)=(Ib-Ib_Prior).*WIb;    
%--------------------------------------------------------------------------------
       
%--------------------------------------------------------------------------------
J=[dFdx dFdy dFdA dFdR dFdb];
%--------------------------------------------------------------------------------