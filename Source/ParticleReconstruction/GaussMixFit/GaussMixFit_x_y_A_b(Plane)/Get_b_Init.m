function b_Init=Get_b_Init(x, y, Ib)
ParticleNum=length(x);
[Ly Lx]=size(Ib);
b_Init=zeros(ParticleNum, 1);
for k=1:ParticleNum      
    yk=y(k); yk=round(yk); yk=min(max(yk,1), Ly);
    xk=x(k); xk=round(xk); xk=min(max(xk,1), Lx);
    b_Init(k)=Ib(yk,xk);    
end

end