function A_Weight=Cal_A_Weight_test(z, b, b_std, MATIRF_Param)
    
Depth=MATIRF_Param.Depth;
I0=MATIRF_Param.I0;
AngleNum=MATIRF_Param.AngleNum;
ParticleNum=length(z);   
A_Weight=ones(ParticleNum, AngleNum);
    
for m=1:AngleNum        
   A_Weight(:,m)=I0(m)*exp(-z/Depth(m))+eps;     
   %A_Weight(:,m)=I0(m)*exp(-z/Depth(m))./(b_std+1)+eps;     
end

for k=1:ParticleNum    
    A_Weight(k,:)=A_Weight(k,:)/sum(A_Weight(k,:));
end