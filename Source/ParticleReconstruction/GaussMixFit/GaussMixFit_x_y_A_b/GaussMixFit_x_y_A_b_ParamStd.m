function A_std=GaussMixFit_x_y_A_R_b_ParamStd(A, PixelNum, F, J)

ParticleNum=length(A);

IdxListP=find(A>0.001);
IdxListN=find(A<=0.001);

IdxList=[IdxListP; ParticleNum+IdxListP;  ParticleNum+IdxListP;  ParticleNum+IdxListP ParticleNum+IdxListP];

subJ=J(:, IdxList);

subF=F(IdxList);

SelectParticleNum=length(IdxListP);

NumDegFree=PixelNum-5*SelectParticleNum;
subParam_std=sqrt(diag(inv(subJ'*subJ))*sum(subF.^2)/NumDegFree);

A_std=zeros(1, ParticleNum);

A_std(IdxListP)=subParam_std(2*ParticleNum+1:3*ParticleNum);
A_std(IdxListN)=max(A_std(IdxListP));