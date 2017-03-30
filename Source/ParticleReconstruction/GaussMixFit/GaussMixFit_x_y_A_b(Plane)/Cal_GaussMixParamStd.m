function [A_Std x_Std y_Std]=Cal_GaussMixParamStd(A, ParticleNum, AngleNum, F, J)

Thresh_A=0.001;
Flag_A=(A>Thresh_A);
Flag_xy=sum(A>Thresh_A, 2);

SelectIdx1=find(Flag_A(:)>0);

subIdx_dFdA=SelectIdx1; 
subIdx_dFdP0=SelectIdx1+AngleNum*ParticleNum;
subIdx_dFdP1=SelectIdx1+AngleNum*2*ParticleNum;
subIdx_dFdP2=SelectIdx1+AngleNum*3*ParticleNum;

SelectIdx2=find(Flag_xy>0);
SelectParticleNum=length(SelectIdx2);

subIdx_dFdx=SelectIdx2+AngleNum*4*ParticleNum;
subIdx_dFdy=SelectIdx2+ParticleNum+AngleNum*4*ParticleNum;

subJ=J(:, [subIdx_dFdA subIdx_dFdP0 subIdx_dFdP1 subIdx_dFdP2 subIdx_dFdx subIdx_dFdy]);

DataNum=length(F);
subNumDegFree=DataNum-(AngleNum*4+2)*SelectParticleNum;
subParam_Std=sqrt(diag(inv(subJ'*subJ))*sum(F.^2)/subNumDegFree);

A_Std=zeros(ParticleNum, AngleNum);
x_Std=zeros(ParticleNum, 1);
y_Std=zeros(ParticleNum, 1);

A_Std(SelectIdx1)=subParam_Std(1:length(SelectIdx1));
x_Std(SelectIdx2)=subParam_Std(end-2*SelectParticleNum+1:end-SelectParticleNum);
y_Std(SelectIdx2)=subParam_Std(end-SelectParticleNum+1:end);

for m=1:AngleNum    
    tempMax=max(A_Std(Flag_A(:,m)>0,m));  
    A_Std(Flag_A(:,m)<1,m)=tempMax;            
end

tempMax_x=max(x_Std(Flag_xy>0));  
x_Std(Flag_xy<1)=tempMax_x;         

tempMax_y=max(y_Std(Flag_xy>0));  
y_Std(Flag_xy<1)=tempMax_y;         