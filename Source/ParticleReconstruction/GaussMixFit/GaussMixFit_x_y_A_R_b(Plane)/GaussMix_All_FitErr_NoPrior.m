function [F J Iobj Ib]=GaussMix_FitErr_NoPrior(Param, ParticleNum, I_Data, X, Y)
% Param = [x(1:ParticleNum) y(1:ParticleNum) a(1:ParticleNum) r(1:ParticleNum) beta(1:ParticleNum)]
% each particle has (x, y, a, r, b)

PixelNum=length(X);

sIdx=1; eIdx=ParticleNum;
x=Param(sIdx:eIdx);

sIdx=ParticleNum+1;   eIdx=sIdx+ParticleNum-1;
y=Param(sIdx:eIdx);

sIdx=2*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
a=Param(sIdx:eIdx);

sIdx=3*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
r=Param(sIdx:eIdx);

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
dFda=zeros(PixelNum,ParticleNum);    
dFdr=zeros(PixelNum,ParticleNum);    
dFdb=zeros(PixelNum,ParticleNum);

for k=1:ParticleNum
    
    tempDist=(x(k)-X).^2+(y(k)-Y).^2;
    
    tempEXP=exp(-tempDist/(2*r(k)^2));      
        
    Iobj=Iobj+a(k)*tempEXP;
        
    b_Weight(:,k)=a(k)*tempEXP+eps;
    Ib=Ib+beta(k).*b_Weight(:,k);      
               
    dFdx(:,k)=a(k)*tempEXP.*((X-x(k))/(r(k)^2));
    
    dFdy(:,k)=a(k)*tempEXP.*((Y-y(k))/(r(k)^2));        
    
    dFda(:,k)=tempEXP;
    
    dFdr(:,k)=a(k)*tempEXP.*tempDist/(r(k)^3);
        
end     

tempSum=sum(b_Weight,2);
for n=1:PixelNum    
    b_Weight(n,:)=b_Weight(n,:)/tempSum(n);
end
Ib=Ib./tempSum;

tempVar=2*(Iobj+Ib)+1;
WI=1./sqrt(tempVar);
%WI=1;

for k=1:ParticleNum    
    
    dFdx(:,k)=dFdx(:,k).*WI;           
    dFdy(:,k)=dFdy(:,k).*WI;   
    dFda(:,k)=dFda(:,k).*WI;                
    dFdr(:,k)=dFdr(:,k).*WI;
                               
    dFdb(:,k)=b_Weight(:,k).*WI;                 
end

%--------------------------------------------------------------------------------
F=(Iobj+Ib-I_Data).*WI;
%--------------------------------------------------------------------------------
    
%--------------------------------------------------------------------------------
sCol=1;  eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFdx;
    
sCol=ParticleNum+1;   eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFdy;
    
sCol=2*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFda;
    
sCol=3*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFdr;    

sCol=4*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFdb;    
%--------------------------------------------------------------------------------