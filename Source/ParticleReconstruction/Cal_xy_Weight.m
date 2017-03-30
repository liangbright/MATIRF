function xy_Weight=Cal_xy_Weight(x_Data, y_Data, A_Data, Astd, Ib_Array, Ib_SigmaArray)
    
[ParticleNum, AngleNum]=size(A_Data);   
[Ly Lx]=size(Ib_Array(:,:,1));   
xy_Weight=ones(ParticleNum, AngleNum);  
for m=1:AngleNum
    
    Ib=Ib_Array(:,:,m);
    Ib_Sigma=Ib_SigmaArray(:,:,m);
    
    for k=1:ParticleNum
        yk=y_Data(k,m); yk=round(yk); yk=min(max(yk,1), Ly);
        xk=x_Data(k,m); xk=round(xk); xk=min(max(xk,1), Lx);
        
        Axy=A_Data(k,m);        
        %bxy=Ib(yk,xk);        
        %bVar=2*Ib(yk,xk);                                            
        bVar=Ib_Sigma(yk,xk).^2;
        
        Var_xy=2*Axy+bVar;
        
        xy_Weight(k,m)=(1-exp(-0.5*Axy^2/Var_xy))+eps;
    end               
end

for k=1:ParticleNum
    xy_Weight(k,:)=xy_Weight(k,:)/sum(xy_Weight(k,:));
end