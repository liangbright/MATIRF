function [b Ib]=Cal_b_plane(x, y, A, R, P0, P1, P2, XPixelList, YPixelList)

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
    
    beta_k=P0(k)+P1(k)*XPixelList+P2(k)*YPixelList;    
    
    Ib=Ib+beta_k.*b_Weight(:,k);                                 
end     

tempSum=sum(b_Weight,2);
Ib=Ib./tempSum;   

b=Ib(IdxList);