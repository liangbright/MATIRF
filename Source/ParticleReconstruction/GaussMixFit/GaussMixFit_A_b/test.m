ParticleNum=10;
x=10+80*rand(ParticleNum,1);
y=10+80*rand(ParticleNum,1);

PSFSigma=1.2;

Ly=100;
Lx=100;
A=300*ones(ParticleNum, 1);
Iobj=GaussMix(PSFSigma, A, x, y, Ly, Lx);
Ib=40;    
I=Iobj+(Iobj*0.02).*randn(Ly,Lx);     
I=poissrnd(I+Ib);
I=I+sqrt(I).*randn(Ly,Lx); 
I=round(I);
%%    
[A_hat_1 b_hat_1]=GaussMix_Fit_A_b_NoPrior_TEST(PSFSigma, x, y, I);

[A_hat b_hat]=GaussMix_Fit_A_b_NoPrior(PSFSigma, x, y, I);
%%
mean(abs(A-A_hat_1))
mean(abs(A-A_hat))

mean(abs(Ib-b_hat_1))
mean(abs(Ib-b_hat))