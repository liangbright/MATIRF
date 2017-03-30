function [F J Ib]=GaussMix_FitErr_NoPrior(Param, ParticleNum, PSFSigma, I_Data, X, Y)
% Param = [A(1:ParticleNum) b(1:ParticleNum) x(1:ParticleNum) y(1:ParticleNum)]
% each particle has (A, b, x, y,)

PixelNum=length(X);

sIdx=1; eIdx=ParticleNum;
A=Param(sIdx:eIdx);

sIdx=ParticleNum+1; eIdx=sIdx+ParticleNum-1;
P0=Param(sIdx:eIdx);

sIdx=2*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
P1=Param(sIdx:eIdx);

sIdx=3*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
P2=Param(sIdx:eIdx);

sIdx=4*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
x=Param(sIdx:eIdx);

sIdx=5*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
y=Param(sIdx:eIdx);
%---------------------------------------------------------------------------------------------
F=zeros(PixelNum,1);
J=zeros(PixelNum,6*ParticleNum);
%%

Iobj=zeros(PixelNum,1);
Ib=zeros(PixelNum,1);
EXPWeight=zeros(PixelNum,ParticleNum);

dFdA=zeros(PixelNum,ParticleNum);    
dFdx=zeros(PixelNum,ParticleNum);    
dFdy=zeros(PixelNum,ParticleNum);    
dFdP0=zeros(PixelNum,ParticleNum);
dFdP1=zeros(PixelNum,ParticleNum);
dFdP2=zeros(PixelNum,ParticleNum);        

for k=1:ParticleNum
    
    tempEXP=exp(-((x(k)-X).^2+(y(k)-Y).^2)/(2*PSFSigma^2));      
        
    Iobj=Iobj+A(k)*tempEXP;
    
    b_k=P0(k)+P1(k)*X+P2(k)*Y;    
    EXPWeight(:,k)=A(k)*tempEXP+eps;
    Ib=Ib+b_k.*EXPWeight(:,k);      
       
    dFdA(:,k)=tempEXP;
        
    dFdx(:,k)=A(k)*((X-x(k))/(PSFSigma^2)).*tempEXP;
    
    dFdy(:,k)=A(k)*((Y-y(k))/(PSFSigma^2)).*tempEXP;
end     

tempSum=sum(EXPWeight,2);
for n=1:PixelNum    
    EXPWeight(n,:)=EXPWeight(n,:)/tempSum(n);
end
Ib=Ib./tempSum;
    
tempVar=2*(Iobj+Ib);
%WI=1./sqrt(tempVar);
WI=1;

for k=1:ParticleNum    
        
    dFdA(:,k)=dFdA(:,k).*WI;            
    dFdx(:,k)=dFdx(:,k).*WI;           
    dFdy(:,k)=dFdy(:,k).*WI;   
                 
    tempProd_WI=EXPWeight(:,k).*WI;              
    dFdP0(:,k)=tempProd_WI;           
    dFdP1(:,k)=X.*tempProd_WI;            
    dFdP2(:,k)=Y.*tempProd_WI;  
end

%--------------------------------------------------------------------------------
F=(Iobj+Ib-I_Data).*WI;
%--------------------------------------------------------------------------------
    
%--------------------------------------------------------------------------------
sCol=1;  eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFdA;
    
sCol=ParticleNum+1;   eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFdP0;
    
sCol=2*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFdP1;
    
sCol=3*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFdP2;    

sCol=4*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFdx;
    
sCol=5*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(:, sCol:eCol)=dFdy;
%--------------------------------------------------------------------------------