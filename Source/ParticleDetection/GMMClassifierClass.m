classdef GMMClassifierClass < handle
    %FEATURECLASSIFIERCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
properties (GetAccess='public',SetAccess='private')

EigenVector

EigenValue

SNR_Thresh

BackGround

NormalizationParam

Type

end
    
methods

%=======================================================================================================================    
function this=GMMClassifierClass()
    
end
%=======================================================================================================================


%=======================================================================================================================
function Label=Build(this, Feature, SNR_Thresh_Init, Type)
    
[Dimension SampleNum]=size(Feature);

[NormalizedFeature Feature_Mean Feature_Sigma]=this.PreNormalization(Feature);

this.NormalizationParam.Feature_Mean=Feature_Mean;
this.NormalizationParam.Feature_Sigma=Feature_Sigma;

if Dimension>1
    [EigenVector, Proj, EigenValue] = Cal_PCA(NormalizedFeature');
    if EigenVector(1,1)<0
        EigenVector=-EigenVector;
        Proj=-Proj;
    end
    
    this.EigenVector=EigenVector;
    this.EigenValue=EigenValue;
       
    ProjFeature=Proj(:,1)';
else
    this.EigenVector=1;
    this.EigenValue=1;
       
    ProjFeature=NormalizedFeature;
end

if nargin == 3
     [Model Feature_SNR]=FitGMM_Auto(this, ProjFeature, SNR_Thresh_Init);
elseif nargin == 4
    if Type == 1    
        [Model Feature_SNR]=FitGMM_Type1(this, ProjFeature, SNR_Thresh_Init);
    else
        [Model Feature_SNR]=FitGMM_Type2(this, ProjFeature, SNR_Thresh_Init);
    end                    
end

Label=Feature_SNR>Model.SNR_Thresh;

this.SNR_Thresh=Model.SNR_Thresh;

this.BackGround.mu=Model.mu;
this.BackGround.Sigma=Model.Sigma;

this.Type=Model.Type;

end
%=======================================================================================================================


%=======================================================================================================================
function [NormalizedFeature Feature_Mean Feature_Sigma]=PreNormalization(this, Feature)

[Dimension SampleNum]=size(Feature);

SNR_Thresh_Init=2;

Feature_Mean=zeros(Dimension, 1);
Feature_Sigma=zeros(Dimension, 1);
Model=cell(Dimension, 1);
for k=1:Dimension
    Model{k}=FitGMM_Type1(this, Feature(k,:), SNR_Thresh_Init);
    Feature_Mean(k)=Model{k}.mu;
    Feature_Sigma(k)=Model{k}.Sigma;
end

NormalizedFeature=bsxfun(@minus, Feature, Feature_Mean);

NormalizedFeature=bsxfun(@rdivide, NormalizedFeature, Feature_Sigma+eps);

end
%=======================================================================================================================


%=======================================================================================================================
function NormalizedFeature=Normalization(this, Feature)

Feature_Mean=this.NormalizationParam.Feature_Mean;
Feature_Sigma=this.NormalizationParam.Feature_Sigma;

NormalizedFeature=bsxfun(@minus, Feature, Feature_Mean);

NormalizedFeature=bsxfun(@rdivide, NormalizedFeature, Feature_Sigma+eps);
end
%=======================================================================================================================


%=======================================================================================================================
function [Model Feature_SNR]=FitGMM_Type1(this, Feature, SNR_Thresh_Init)
    
Fit_Init= gmdistribution.fit(Feature', 1);
mu_Init=Fit_Init.mu;
Sigma_Init=sqrt(Fit_Init.Sigma(:,:,1));

mu=mu_Init;
Sigma=Sigma_Init;

tempFeature=Feature;

OK=0;
counter=0;
while OK==0 && counter<3
    
    counter=counter+1;
    
    tempFeature=tempFeature(tempFeature<mu+2.5*Sigma);

    Fit=gmdistribution.fit(tempFeature', 1);

    if abs(mu-Fit.mu) < Sigma/10
        OK=1;
    else    
        mu=Fit.mu;
        Sigma=sqrt(Fit.Sigma(:,:,1));
    end
end

Model.mu=mu;
Model.Sigma=Sigma;
Model.SNR_Thresh=SNR_Thresh_Init;
Model.Type=1;

Feature_SNR=(Feature-mu)/Sigma;
%pValue=1-normcdf(Feature, mu, Sigma);    

end 
%=======================================================================================================================


%=======================================================================================================================
function [Model Feature_SNR]=FitGMM_Type2(this, Feature, SNR_Thresh_Init)

Model=[];
Feature_SNR=[];

Fit=gmdistribution.fit(Feature',2);

if Fit.Converged >0
    mu_a=Fit.mu(1);
    Sigma_a=sqrt(Fit.Sigma(:,:,1));
    A_a= Fit.PComponents(1);

    mu_b=Fit.mu(2);
    Sigma_b=sqrt(Fit.Sigma(:,:,2));
    A_b= Fit.PComponents(2);
   
    Thresh=this.Cal_Thresh(mu_a, Sigma_a, mu_b, Sigma_b, A_a, A_b);
    
    if mu_a < mu_b
        mu=mu_a;
        Sigma=Sigma_a;
    else
        mu=mu_b;
        Sigma=Sigma_b;
    end
        
    if mu < Thresh
        SNR=(Thresh-mu)/Sigma;        
    else
        SNR=0;
    end
    
    if SNR > SNR_Thresh_Init        
        Model.mu=mu;
        Model.Sigma=Sigma;
        Model.SNR_Thresh=SNR;      
        Model.Type=2;
        Feature_SNR=(Feature-mu)/Sigma;
        %Feature_pValue=1-normcdf(Feature, mu, Sigma);    
    %else
    % Model=[];    
    end   
%else
%   Model=[];
end

end
%=======================================================================================================================  
    

%=======================================================================================================================
function [Model Feature_SNR]=FitGMM_Auto(this, Feature, SNR_Thresh_Init)
    
[Model_1 Feature_SNR_1]=this.FitGMM_Type1(Feature, SNR_Thresh_Init);

[Model_2 Feature_SNR_2]=this.FitGMM_Type2(Feature, SNR_Thresh_Init);

if isempty(Model_2)    
    %ClusterNum=1;
    Model=Model_1;
    Feature_SNR=Feature_SNR_1;
    disp('1 Class')
else
    %ClusterNum=2;
    Model=Model_2;
    Feature_SNR=Feature_SNR_2;
    disp('2 Classes')
end

end 
%=======================================================================================================================


%=======================================================================================================================
function Thresh=Cal_Thresh(this, m1, s1, m2, s2, p1, p2)

% ax^2+bx+c=0

R=(s2/s1)^2;
a=0.5*(R-1);
b=-R*m1+m2;
c=0.5*(R*m1^2-m2^2)-s2^2*log(p1*s2/(p2*s1));

temp=sqrt(b^2-4*a*c);
if ~isreal(temp)
    Thresh=[];
end

Thresh1=(-b+temp)/(2*a);
Thresh2=(-b-temp)/(2*a);

max_m=max(m1,m2);
min_m=min(m1,m2);
if Thresh1 >=min_m && Thresh1 <= max_m
    Thresh=Thresh1;
elseif Thresh2 >=min_m && Thresh2 <= max_m
    Thresh=Thresh2;
else
    Thresh=[];
end

end
%=======================================================================================================================


%=======================================================================================================================
function [Label SNR pValue]=Identify(this, Feature)

NormalizedFeature=this.Normalization(Feature);

Proj=NormalizedFeature'*this.EigenVector;
ProjFeature=Proj(:,1)';

mu=this.BackGround.mu;
Sigma=this.BackGround.Sigma;

SNR=(ProjFeature-mu)/Sigma;
Label=SNR>this.SNR_Thresh;

%pValue=1-normcdf(ProjFeature, mu, Sigma);    
pValue=1-normcdf(ProjFeature, mu, Sigma);    

end
%=======================================================================================================================

end
end