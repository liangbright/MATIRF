function R_Weight=Cal_R_Weight(x, y, A_Data, A_DataSigma, Ib_Array, Ib_SigmaArray)
    
[ParticleNum, AngleNum]=size(A_Data);    
[Ly Lx]=size(Ib_Array(:,:,1));
R_Weight=ones(ParticleNum, AngleNum);    
for m=1:AngleNum
    
    Ib=Ib_Array(:,:,m);
    Ib_Sigma=Ib_SigmaArray(:,:,m);
    
    for k=1:ParticleNum
        
        yk=y(k); yk=round(yk); yk=min(max(yk,1), Ly);
        xk=x(k); xk=round(xk); xk=min(max(xk,1), Lx);
    
        Axy=A_Data(k,m);
        %Asigma=A_DataSigma(k,m);        
        
        %bxy=Ib(yk,xk);         
        %bVar=2*Ib(yk,xk);
        bVar=Ib_Sigma(yk,xk)^2;
        
        Var_xy=2*Axy+bVar;                
        
        R_Weight(k,m)=(1-exp(-0.5*Axy^2/Var_xy))+eps;             
    end    
end

for k=1:ParticleNum    
    R_Weight(k,:)=R_Weight(k,:)/sum(R_Weight(k,:));
end