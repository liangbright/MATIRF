function [b Ib]=Cal_b(x, y, A, R, beta, XPixelList, YPixelList)

ParticleNum=length(x);
PixelNum=length(XPixelList);

IdxList=zeros(ParticleNum, 1);

Ib=zeros(PixelNum,1);
b_Weight=zeros(PixelNum,ParticleNum);

for k=1:ParticleNum
    
    tempDist=(x(k)-XPixelList).^2+(y(k)-YPixelList).^2;
    
    [~, tempIndex]=min(tempDist);
    
    IdxList(k)=tempIndex;
    
    tempEXP=exp(-tempDist/(2*R(k)^2));                  
        
    b_Weight(:,k)=A(k)*tempEXP+eps;
    Ib=Ib+beta(k).*b_Weight(:,k);                                 
end     

tempSum=sum(b_Weight,2);
Ib=Ib./tempSum;   

b=Ib(IdxList);