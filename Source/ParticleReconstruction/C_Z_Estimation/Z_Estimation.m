function z=Z_Estimation(C, Range_z, A_Data, A_Weight, MATIRF_Param)

[ParticleNum, AngleNum]=size(A_Data);

z=zeros(ParticleNum,1);

I0=MATIRF_Param.I0;
Depth=MATIRF_Param.Depth;
%%
Options = optimset('Jacobian','on','Display','off');

Flag=zeros(ParticleNum,1);

Amp_Thresh=1;

for k=1:ParticleNum
                    
    z_min=Range_z(1, k);
    z_max=Range_z(2, k);
        
    tempAk=(A_Data(k,:))';
        
    Index=find(tempAk>Amp_Thresh);
    if length(Index) >= 1
            
        z_init=z_min;                    
        Lowbd=z_min;
        Upbd=z_max;    
        W=A_Weight(k,:)';                    
        z(k)= lsqnonlin(@DecayProfile_FitErr_z, z_init, Lowbd, Upbd, Options, C(k), tempAk, I0, Depth, W);
    else
        z(k)=z_max; 
        Flag(k)=1;        
    end
end

z_mean=mean(z(Flag<1));
if ~isempty(z_mean)
    z(Flag==1)=z_mean;
end
    