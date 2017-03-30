PSFSigma=1.2;
Ly=110;
Lx=20;
[X Y]=meshgrid(1:Lx, 1:Ly);

y0=(10:10:100)'+0.5;
x0=10*ones(10,1)+0.5;
ParticleNum=length(x0);

z0=ones(ParticleNum,1);
R0=(25:(25/(ParticleNum-1)):50)';

Ibg=40*ones(Ly,Lx);

Iobj=GaussMix(PSFSigma, 300*ones(ParticleNum,1), x0, y0, Lx, Ly);

I=Iobj+(Iobj*0.02).*randn(Ly,Lx);     
I=poissrnd(I+Ibg);
I=I+sqrt(I).*randn(Ly,Lx); 
I=round(I);
%%
x_Init=round(x0);
y_Init=round(y0);

tempA=GaussMix_Fit_A_b(PSFSigma, x_Init, y_Init, I);        
Iobj_Init=GaussMix(PSFSigma, tempA, x_Init, y_Init, Lx, Ly);                
I_Weight=(Iobj_Init<0.1);      
[Ib Ib_Sigma]=BackgroundEstimation(I, I_Weight, PSFSigma, 'B');        
A_Init=GaussMix_Fit_A_b(PSFSigma, x_Init, y_Init, I, Ib, Ib_Sigma);            
        
%%    
b_Init=zeros(ParticleNum, 1);
for k=1:ParticleNum      
    yk=y_Init(k); yk=round(yk); yk=min(max(yk,1), Ly);
    xk=x_Init(k); xk=round(xk); xk=min(max(xk,1), Lx);
    b_Init(k)=Ib(yk,xk);    
end
%% Fit    
[Param A x y b Ic RegionInfo Region F J] = GaussMix_Fit(PSFSigma, A_Init, x_Init, y_Init, b_Init, I);
%%
[Param A x y b Ic RegionInfo Region F J] = GaussMix_Fit(PSFSigma, A_Init, x_Init, y_Init, b_Init, I, Ib, Ib_Sigma);
%%
[A_Std x_Std y_Std]=Cal_GaussMixParamStd(A, ParticleNum, 1, F, J);
%%
[A_hat b_hat A_Std_hat b_Std_hat]=GaussMix_Fit_A_b(PSFSigma, x, y, I);  
%%
%[A_hat_2 b_hat_2 A_Std_hat_2 b_Std_hat_2]=GaussMix_Fit_A_b(PSFSigma, x, y, I, Ib, Ib_Sigma, A, 0.001*A);

%%
A_RMSE=sqrt(sum((A-300).^2)/ParticleNum);
x_RMSE=sqrt(sum((x-x0).^2)/ParticleNum);
y_RMSE=sqrt(sum((y-y0).^2)/ParticleNum);
b_RMSE=sqrt(sum((b-40).^2)/ParticleNum);
