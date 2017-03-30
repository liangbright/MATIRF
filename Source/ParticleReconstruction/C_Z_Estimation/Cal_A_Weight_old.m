function A_Weight=Cal_A_Weight_old(x, y, A_Data, A_DataSigma, Ib_Array, Ib_SigmaArray)
    
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
        
        ISigma=sqrt(2*(Asigma+bSigma/AngleNum));
         
        A_Weight(k,m)=1./ISigma;
        
        %A_Weight(k,m)=1;
        
        %A_Weight(k,m)=A/sqrt(2*(A+b));
        
        %A_Weight(k,m)=1-exp(-A^2/(2*b));
        
        %A_Weight(k,m)=1-exp(-A^2/(2*(A+b)));
        
        %A_Weight(k,m)=exp(-0.5*(A/Asigma^2));
        
        %A_Weight(k,m)=(1-exp(-A^2/(2*bxy)))/Asigma;
                        
        %A_Weight(k,m)=(1-exp(-A^2/bVar))/Asigma;  
        
        %A_Weight(k,m)=1-exp(-(A/ISigma)^2);  
        
        %A_Weight(k,m)=1/(Asigma+1);
        
        %A_Weight(k,m)=1-exp(-0.5*A^2/Var_xy))+eps;  
                
    end    
end

%A_Weight=sqrt(A_Weight);
%{
for k=1:ParticleNum    
    A_Weight(k,:)=A_Weight(k,:)/sum(A_Weight(k,:));
end
%}

A_Weight=A_Weight/max(A_Weight(:));