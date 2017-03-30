function A_Weight=Cal_A_Weight(x, y, A_Data, A_DataSigma, Ib_Array, Ib_SigmaArray)
    
[ParticleNum, AngleNum]=size(A_Data);
    
[Ly Lx]=size(Ib_Array(:,:,1));

A_Weight=ones(ParticleNum, AngleNum);

for m=1:AngleNum
    
    Ib=Ib_Array(:,:,m);
    Ib_Sigma=Ib_SigmaArray(:,:,m);
    
    for k=1:ParticleNum
        
        yk=y(k); yk=round(yk); yk=min(max(yk,1), Ly);
        
        xk=x(k); xk=round(xk); xk=min(max(xk,1), Lx);
    
        A=A_Data(k,m);
        Asigma=A_DataSigma(k,m);        
        
        b=Ib(yk,xk);         
        bSigma=Ib_Sigma(yk,xk);
        
        ISigma=sqrt(2*(A+b));
        
        %ISigma=Asigma+bSigma;
        
        A_Weight(k,m)=1./ISigma;
                
    end    
end

A_Weight=A_Weight/max(A_Weight(:));