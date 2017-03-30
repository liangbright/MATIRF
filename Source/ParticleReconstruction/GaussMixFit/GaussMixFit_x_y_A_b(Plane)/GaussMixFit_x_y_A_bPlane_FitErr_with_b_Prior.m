function [F J Ib]=GaussMix_FitErr_with_b_Prior(Param, ParticleNum, PSFSigma, I_Data, X, Y, ...
                                               Ib_Prior, Ib_SigmaPrior)
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
F=zeros(2*PixelNum,1);
J=zeros(2*PixelNum,6*ParticleNum);
%%

Iobj=zeros(PixelNum,1);
Ib=zeros(PixelNum,1);
EXPWeight=zeros(PixelNum,ParticleNum);

dFdA=zeros(2*PixelNum,ParticleNum);    
dFdx=zeros(2*PixelNum,ParticleNum);    
dFdy=zeros(2*PixelNum,ParticleNum);    
dFdP0=zeros(2*PixelNum,ParticleNum);
dFdP1=zeros(2*PixelNum,ParticleNum);
dFdP2=zeros(2*PixelNum,ParticleNum);        

for k=1:ParticleNum
    
    tempEXP=exp(-((x(k)-X).^2+(y(k)-Y).^2)/(2*PSFSigma^2));      
        
    Iobj=Iobj+A(k)*tempEXP;
    
    b_k=P0(k)+P1(k)*X+P2(k)*Y;    
    EXPWeight(:,k)=A(k)*tempEXP+eps;
    Ib=Ib+b_k.*EXPWeight(:,k);      
       
    dFdA(1:PixelNum,k)=tempEXP;
        
    dFdx(1:PixelNum,k)=A(k)*((X-x(k))/(PSFSigma^2)).*tempEXP;
    
    dFdy(1:PixelNum,k)=A(k)*((Y-y(k))/(PSFSigma^2)).*tempEXP;
end     

tempSum=sum(EXPWeight,2);
for n=1:PixelNum    
    EXPWeight(n,:)=EXPWeight(n,:)/tempSum(n);
end
Ib=Ib./tempSum;

tempVar=2*(Iobj+Ib);
WI=1./sqrt(tempVar); 

WIb=1./Ib_SigmaPrior;
    
for k=1:ParticleNum    
        
    dFdA(1:PixelNum,k)=dFdA(1:PixelNum,k).*WI;            
    dFdx(1:PixelNum,k)=dFdx(1:PixelNum,k).*WI;           
    dFdy(1:PixelNum,k)=dFdy(1:PixelNum,k).*WI;   
                 
    tempProd_WI=EXPWeight(:,k).*WI;              
    dFdP0(1:PixelNum,k)=tempProd_WI;           
    dFdP1(1:PixelNum,k)=X.*tempProd_WI;            
    dFdP2(1:PixelNum,k)=Y.*tempProd_WI;         
    
    % Cost=sum{(Ib-b_Data).^2}
    tempProd_WIb=EXPWeight(:,k).*WIb;    
    sIdx=PixelNum+1; eIdx=2*PixelNum;
    dFdP0(sIdx:eIdx,k)=tempProd_WIb;    
    dFdP1(sIdx:eIdx,k)=X.*tempProd_WIb;   
    dFdP2(sIdx:eIdx,k)=Y.*tempProd_WIb;   
end
    
%--------------------------------------------------------------------------------
F(1:PixelNum)=(Iobj+Ib-I_Data).*WI;
F(PixelNum+1:2*PixelNum)=(Ib-Ib_Prior).*WIb;    
%--------------------------------------------------------------------------------
    
%--------------------------------------------------------------------------------
sRow=1; eRow=2*PixelNum;    
%--------------------------------------------------------------------------------
sCol=1;  eCol=sCol+ParticleNum-1;
J(sRow:eRow, sCol:eCol)=dFdA;
    
sCol=ParticleNum+1;   eCol=sCol+ParticleNum-1;
J(sRow:eRow, sCol:eCol)=dFdP0;
    
sCol=2*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(sRow:eRow, sCol:eCol)=dFdP1;
    
sCol=3*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(sRow:eRow, sCol:eCol)=dFdP2;    

sCol=4*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(sRow:eRow, sCol:eCol)=dFdx;
    
sCol=5*ParticleNum+1; eCol=sCol+ParticleNum-1;
J(sRow:eRow, sCol:eCol)=dFdy;
%--------------------------------------------------------------------------------