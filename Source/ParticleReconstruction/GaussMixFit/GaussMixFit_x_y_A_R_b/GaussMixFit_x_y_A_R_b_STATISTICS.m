function [A_std MSE]=GaussMixFit_x_y_A_R_b_STATISTICS(Param, ParticleNum, I_Data, XPixelList, YPixelList)

[F J Iobj_model Ib_model]=GaussMixFit_x_y_A_R_b_FitErr_NoPrior(Param, ParticleNum, I_Data, XPixelList, YPixelList);
%%
MSE=mean((Iobj_model+Ib_model-I_Data).^2); 
%%

sIdx=2*ParticleNum+1; eIdx=sIdx+ParticleNum-1;
A=Param(sIdx:eIdx);
A=A(:)';

A_std=zeros(ParticleNum, 1);

IdxListP=find(A>0.001);
IdxListN=find(A<=0.001);

if ~isempty(IdxListP)
    
    IdxList=[IdxListP ParticleNum+IdxListP 2*ParticleNum+IdxListP 3*ParticleNum+IdxListP 4*ParticleNum+IdxListP];
    subJ=J(:, IdxList);

    PixelNum=length(XPixelList);
    subF=F(1:PixelNum);

    SelectParticleNum=length(IdxListP);

    NumDegFree=PixelNum-5*SelectParticleNum;
    subParam_std=sqrt(diag(inv(subJ'*subJ))*sum(subF.^2)/NumDegFree);

    A_std(IdxListP)=subParam_std(2*ParticleNum+1:3*ParticleNum);
    A_std(IdxListN)=max(A_std(IdxListP));
else
    A_std(:)=NaN; 
end