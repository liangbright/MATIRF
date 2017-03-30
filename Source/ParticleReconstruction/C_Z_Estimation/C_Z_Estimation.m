function [C z A A_Sigma] = C_Z_Estimation(Range_C, Range_z, A_Data, A_DataSigma, A_Weight, MATIRF_Param, LoopIndex)

[ParticleNum, AngleNum]=size(A_Data);

C=zeros(ParticleNum,1);
z=zeros(ParticleNum,1);
A=zeros(ParticleNum,AngleNum);
%A_Sigma=zeros(ParticleNum,AngleNum);

I0=MATIRF_Param.I0;
Depth=MATIRF_Param.Depth;
%MAPE_Depth=MATIRF_Param.MAPE_Depth;

%%
Options = optimset('Jacobian','on','Display','off');

Flag=zeros(ParticleNum,1);

Amp_Thresh=1;

for k=1:ParticleNum
    
    tempAk=A_Data(k,:)';       
    
    Adiff=diff([-1; tempAk]);
    IdxList=find(Adiff>0);
    
    if length(IdxList) < AngleNum
        disp(['weired A data in  C_Z_Estimation: Particle ID ' num2str(k), ', A ' num2str(tempAk')]);
    end

    IdxList_inti=IdxList(tempAk(IdxList)>Amp_Thresh);
     
    if length(IdxList_inti) >= 2
    
        [C_init z_init]=Estimate_InitParam(tempAk(IdxList_inti), I0(IdxList_inti), Depth(IdxList_inti));
                        
        C_init=max(C_init, 1);
        z_init=max(z_init, 1);      
        
        z_min=Range_z(1, k);
        z_max=Range_z(2, k);
        
        C_min=Range_C(1, k);
        C_max=Range_C(2, k);
        
        Init= [C_init; z_init];
        Lowbd=[C_min;  z_min];
        Upbd= [C_max;  z_max];
    
        W=A_Weight(k,:)';        
         
       % IdxList_fit=find(tempAk>Amp_Thresh);
                  
       % IdxList_fit=IdxList_inti;
        
        IdxList_fit=1:length(tempAk);
       
        Param=lsqnonlin(@DecayProfile_FitErr, Init, Lowbd, Upbd, Options, ...
                        tempAk(IdxList_fit), I0(IdxList_fit), Depth(IdxList_fit), W(IdxList_fit));
              
        C(k)=Param(1);
        z(k)=Param(2);     
        
    else
        Flag(k)=1;
    end       
    
end

C_mean=mean(C(Flag<1));
z_mean=mean(z(Flag<1));

for k=1:ParticleNum
            
    if Flag(k)>0
        
        z_min=Range_z(1, k);
        z_max=Range_z(2, k);
        
        tempAk=(A_Data(k,:))';
        
        Index=find(tempAk>Amp_Thresh);
        if length(Index) >= 1
            C(k)=C_mean;           
            
            z_init=z_min;                    
            Lowbd=z_min;
            Upbd=z_max;    
            W=A_Weight(k,:)';         

            z(k)= lsqnonlin(@DecayProfile_FitErr_z, z_init, Lowbd, Upbd, Options, C(k), tempAk, I0, Depth, W);
        else % 0
            C(k)=C_mean;
            %z(k)=z_mean;
            z(k)=z_max;
        end
    end
end
    
%%
for k=1:ParticleNum
    A(k,:)=C(k).*DecayProfile(z(k), I0, Depth)';
end
%%
%{
NumDegFree=AngleNum-2;
NumDegFree=max(NumDegFree, 1);

A_FitSigma=ones(ParticleNum, AngleNum);
%A_ModelSigma=ones(ParticleNum, AngleNum);
for k=1:ParticleNum
    
    Residual=(A(k,:)-A_Data(k,:)).*A_Weight(k,:);
    
    A_FitSigma(k,:)=sqrt(sum(Residual.^2)/NumDegFree);
    
   % A_ModelSigma(k,:)=(MAPE_Depth*0.00375)*I0*C(k);
end

%A_Sigma=sqrt(0*A_DataSigma.^2+0*A_ModelSigma.^2+A_FitSigma.^2);
%A_Sigma=sqrt((A_DataSigma.*A_FitSigma).^2+A_ModelSigma.^2);
%A_Sigma=A_DataSigma.*A_FitSigma;
%A_Sigma=0.02+A_DataSigma.*min(A_FitSigma, 1/Loop);
%A_Sigma=0.02+A_DataSigma/Loop;
%A_Sigma=A_ModelSigma+A_DataSigma*(1-exp(-1/Loop));
%A_Sigma=A_ModelSigma+A_DataSigma/Loop;
%}
A_Sigma=0.01+A_DataSigma/sqrt(LoopIndex);