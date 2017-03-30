Depth=[100 120 150 200 250 300 350 400 450 inf]';
I0=[1 1 1 1 1 1 1 1 1 0.1]';
%%
Depth=[10 100]';
I0=[1 1 ]';
C=100;
z=500;

W=ones(length(I0),1);
%%
for k=1:100
A=C*exp(-z./Depth); %+(C*0.01).*randn(size(I0)); 
A(1)=0;
%%
Init= [0; 0];
Lowbd=[10; 10];
Upbd= [1000; 1000];
Options = optimset('Jacobian','on','Display','off');

Param=Init;
Param= lsqnonlin(@DecayProfile_FitErr, Param, Lowbd, Upbd, Options, ...
                  A, I0, Depth, W);
%%
A_Data=A';
[C_hat z_hat] = Module_2(A_Data, I0, Depth, W);
%%              
%C_hat=Param(1);
%z_hat=Param(2);     

C_err(k)=C_hat-Signal;
z_err(k)=z_hat-z;

end

 k=1:100;
figure(1); plot(k, C_err, 'bo-', k, z_err, 'ro-'); grid on
