function Iobj=GaussMix(x, y, A, R, Lx, Ly, GaussMixtureRange_Thresh)

if nargin < 7
    Amp_Thresh_log=-6.2146; %log(0.002);
else
    Amp_Thresh_log=log(min(0.002, GaussMixtureRange_Thresh));
end
% fast  1s/100
%--------------------------------------------
Iobj=zeros(Ly, Lx);
ParticleNum=length(A);

for k=1:ParticleNum    
        
    tempA_k=max(A(k), 0.002);
    delta=1.4142*R(k)*(log(tempA_k)-Amp_Thresh_log);    
    delta=ceil(delta);
    
    x1=round(x(k)-delta); x1=max(x1, 1);
    x2=round(x(k)+delta); x2=min(x2, Lx);
    y1=round(y(k)-delta); y1=max(y1, 1);
    y3=round(y(k)+delta); y3=min(y3, Ly);
    
    [X Y]=meshgrid(x1:x2, y1:y3);                    
    tempI=A(k)*exp(-((x(k)-X).^2+(y(k)-Y).^2)/(2*R(k)^2));    
    Iobj(y1:y3, x1:x2)=Iobj(y1:y3, x1:x2)+tempI;     
end

%{
[X Y]=meshgrid(1:Lx, 1:Ly);
for k=1:ParticleNum                   
    tempEXP=exp(-((x(k)-X).^2+(y(k)-Y).^2)/(2*R(k)^2));    
    Iobj=Iobj+A(k)*tempEXP;     
end
%}