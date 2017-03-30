function [A_Init b_Init Radius_xy Range_A Range_R Range_b]=GaussMixFit_InitialParam(x_Init, y_Init, RaidusMaxRange, PSFSigma, I, ObjLabel, Ib, Ib_Sigma)

ParticleNum=length(x_Init);

A_Init=zeros(1, ParticleNum);
b_Init=zeros(1, ParticleNum);
Radius_xy=zeros(1, ParticleNum);
Range_A=zeros(2, ParticleNum);
Range_R=zeros(2, ParticleNum);
Range_b=zeros(2, ParticleNum);

[Ly Lx]=size(I);

Imin=min(I(ObjLabel>0));
    
%-------------------------------------------------------------------------------
Radius_xy(:)=PSFSigma+0.5;

%-------------------------------------------------------------------------------
if ~isempty(Ib)
    for k=1:ParticleNum
        yk=y_Init(k); yk=round(yk); yk=min(max(yk,1), Ly);
        xk=x_Init(k); xk=round(xk); xk=min(max(xk,1), Lx);
    
        b_Init(k)=Ib(yk,xk);    
        A_Init(k)=max([I(yk,xk)-b_Init(k), Ib_Sigma(yk,xk)]);
    end
else
    for k=1:ParticleNum
        yk=y_Init(k); yk=round(yk); yk=min(max(yk,1), Ly);
        xk=x_Init(k); xk=round(xk); xk=min(max(xk,1), Lx);
       
        A_Init(k)=max([I(yk,xk)-Imin, 1]);
    end
    
    b_Init(:)=Imin;
end

Range_A(1,:)=0;
Range_A(2,:)=5*A_Init+median(A_Init);

%-------------------------------------------------------------------------------
Range_R(1,:)=RaidusMaxRange(1);
Range_R(2,:)=RaidusMaxRange(2);

%-------------------------------------------------------------------------------
if ~isempty(Ib)
    for k=1:ParticleNum
        yk=y_Init(k); yk=round(yk); yk=min(max(yk,1), Ly);
        xk=x_Init(k); xk=round(xk); xk=min(max(xk,1), Lx);
        bk=Ib(k);
        bk_std=Ib_Sigma(yk,xk);
        
        Range_b(1, k)=max(bk-3*bk_std-1, 0);
        Range_b(2, k)=bk+3*bk_std+1;
    end
else
    
    Range_b(1,:)=0;
    Range_b(2,:)=Range_A(2,:);

end