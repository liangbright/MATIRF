function [Region AllRegionInfo RegionIdxList RegionLabel]=SelectGaussMixRegion(I, ObjLabel, x, y, A, R, b, MaxWindowRadius)

if nargin < 8
    MaxWindowRadius=9;
end

dist_ratio_Thresh=0.9;

[Ly, Lx]=size(I);

A=max(A, 1);
R=max(R, 1);

ParticleNum=length(A);

Iobj=zeros(Ly, Lx);
tempRegionLabel=zeros(Ly, Lx);
tempLabelNum=0;

tempRegionIdxList=zeros(ParticleNum, 1);

%b_Sigma=sqrt(2*b);
%Amp_Thresh=max(b_Sigma, 0.002);    
Amp_Thresh=0.002;
alpha=1.4142*sqrt(log(A)-log(Amp_Thresh));
alpha=max(alpha, 3);
alpha=min(alpha, 4); 
    
delta=alpha.*R;
delta=min(delta, MaxWindowRadius);
delta=ceil(delta);
    
%test the effect of the window size
%delta=6*ones(ParticleNum, 1);

for k=1:ParticleNum
    
    x1=round(x(k)-delta(k)); x1=max(x1, 1);
    x2=round(x(k)+delta(k)); x2=min(x2, Lx);
    y1=round(y(k)-delta(k)); y1=max(y1, 1);
    y3=round(y(k)+delta(k)); y3=min(y3, Ly);
    
    [X Y]=meshgrid(x1:x2, y1:y3);      
    dist=sqrt((x(k)-X).^2+(y(k)-Y).^2);
    tempI=A(k)*exp(-dist.^2/(2*R(k)^2));    
    
    tempIobj=Iobj(y1:y3, x1:x2);
    
    IdxList=find(tempIobj>0 & dist < delta(k)*dist_ratio_Thresh);   
    
    if isempty(IdxList)                
        tempIdxList=find(tempI>tempIobj & dist<=delta(k));
        X_k=X(tempIdxList);
        Y_k=Y(tempIdxList);        
        IdxList_k= sub2ind([Ly Lx], Y_k, X_k);        
        tempLabelNum=tempLabelNum+1;
        tempRegionLabel(IdxList_k)=tempLabelNum;      
        tempRegionIdxList(k)=tempLabelNum;
    else
        tempX=X(IdxList);
        tempY=Y(IdxList);     
        tempIdxList= sub2ind([Ly Lx], tempY, tempX);        
        Idx_Label=tempRegionLabel(tempIdxList);
        Idx_Label=unique(Idx_Label(Idx_Label>0));
        
        X_k=X(dist<=delta(k));
        Y_k=Y(dist<=delta(k));
        IdxList_k= sub2ind([Ly Lx], Y_k, X_k);        
        
        tempLabelNum=tempLabelNum+1;       
        tempRegionLabel(IdxList_k)=tempLabelNum;
        tempRegionIdxList(k)=tempLabelNum;
         
        for n=1:length(Idx_Label)
            tempRegionLabel(tempRegionLabel==Idx_Label(n))=tempLabelNum;
            tempRegionIdxList(tempRegionIdxList==Idx_Label(n))=tempLabelNum;
        end
        
    end                
    
    Iobj(y1:y3, x1:x2)=Iobj(y1:y3, x1:x2)+tempI;                    
end
%%
tempRegionLabel(ObjLabel==0)=0;

tempRegion= regionprops(tempRegionLabel, 'PixelList','PixelIdxList');
tempRegionNum=length(tempRegion);

if tempRegionNum ~= tempLabelNum
    disp('Wrong in SelectGaussMixRegion');
end

RegionIdxList=zeros(ParticleNum, 1);
RegionLabel=zeros(Ly, Lx);
Region=[];
RegionNum=0;
for n=1:tempRegionNum
    if ~isempty(tempRegion(n).PixelIdxList)        
        
        RegionNum=RegionNum+1;
        Region(RegionNum).RegionIdx=RegionNum;
        
        Region(RegionNum).ParticleIdx=[];
        Region(RegionNum).x=[];
        Region(RegionNum).y=[];         
        Region(RegionNum).A=[];
        Region(RegionNum).R=[];
        Region(RegionNum).b=[];    
        
        PixelIdxList=tempRegion(n).PixelIdxList;
        Region(RegionNum).PixelIdxList=PixelIdxList;
        Region(RegionNum).XPixelList=tempRegion(n).PixelList(:,1);
        Region(RegionNum).YPixelList=tempRegion(n).PixelList(:,2);    
        
        RegionLabel(PixelIdxList)=RegionNum;                
        
        RegionIdxList(tempRegionIdxList==n)=RegionNum;
        
    end
end

for k=1:ParticleNum
    
    Index=RegionIdxList(k);       
    
    Region(Index).ParticleIdx(end+1,1)=k;        
    Region(Index).x(end+1,1)=x(k);
    Region(Index).y(end+1,1)=y(k);    
    Region(Index).A(end+1,1)=A(k);
    Region(Index).R(end+1,1)=R(k);
    Region(Index).b(end+1,1)=b(k);                
    
end

AllRegionInfo.PixelIdxList=[];
AllRegionInfo.XPixelList=[];
AllRegionInfo.YPixelList=[];
for m=1:RegionNum
    tempNum=length(Region(m).PixelIdxList);
    if tempNum > 0
        AllRegionInfo.PixelIdxList(end+1:end+tempNum, 1)=Region(m).PixelIdxList;
        AllRegionInfo.XPixelList(end+1:end+tempNum, 1)=Region(m).XPixelList;
        AllRegionInfo.YPixelList(end+1:end+tempNum, 1)=Region(m).YPixelList;
    end
end