function FeaturePoint=FeaturePointDetection(I, CellLabel, ParticleRadius)

[Ly Lx]=size(I);

%Hfilter_LoG=Mask_LoG(ParticleRadius);
%I_LoG=imfilter(I, Hfilter_LoG, 'symmetric');

Hfilter_MF=Mask_MatchedFilter(ParticleRadius);
I_MF=imfilter(I, Hfilter_MF, 'symmetric');

Hfilter_TopHat=Mask_TopHat(ParticleRadius);    
I_TPH=imfilter(I, Hfilter_TopHat, 'symmetric');

%%
%Is=I_LoG;
Is=I_MF;
%% find local maxima and extract feature

maskAsize=2*round(ParticleRadius)+1;
maskAsize=max(3, maskAsize);
maskA=ones(maskAsize, maskAsize);
[maskX,maskY] = meshgrid(1:maskAsize,1:maskAsize);
dist=sqrt((maskX-round(maskAsize/2)).^2+(maskY-round(maskAsize/2)).^2);
maskA(dist>((maskAsize)/2))=0;
numA=sum(maskA(:));

Isp=Is; 
Isp(Isp<0)=0;

Isn=-Is; 
Isn(Isn<0)=0;
    
Iord_p = ordfilt2(Isp,numA,maskA);
Iord_n = ordfilt2(Isn,numA,maskA);
Iord_p=Iord_p.*(Iord_p>Iord_n).*(Isp>0);

In_a=Isp;
In_a(Iord_p ~= Isp) = 0;
In_a(I_TPH<0)=0;

In_b=I_TPH;
%In_b(In_a==0)=0;

%% Output Feature Points
IndexList=find(In_a>0 & CellLabel>0);
PointNum=length(IndexList);

if PointNum >= 1

    FeaturePoint=zeros(4, PointNum);

    [X Y]=meshgrid(1:Lx, 1:Ly);
    FeaturePoint(1,:)=X(IndexList);
    FeaturePoint(2,:)=Y(IndexList);
    FeaturePoint(3,:)=In_a(IndexList);    
    FeaturePoint(4,:)=In_b(IndexList);
    
else
    FeaturePoint=[];
end
