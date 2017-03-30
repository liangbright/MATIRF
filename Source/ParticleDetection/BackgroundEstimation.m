function [b b_Sigma]= BackgroundEstimation(I, SegLabel, PSFSigma, TypicalRadius)

% SegLabel  :   0 - out of the cell, 1 - membrane background, 2 - Gausian Mixture of the particles

% maximum iteration number
MaxIter=8;

% Window Size for local estimation
% for method A
R=5*TypicalRadius;
% for method B
%R=6*TypicalRadius;

R=max(R, 4*PSFSigma);
R=min(R, 6*PSFSigma);
R=round(R);

b_Label=(SegLabel==1);

Obj_Label=(SegLabel>0);

b=LocalEstimation_A(I, b_Label, Obj_Label, R, MaxIter);

%b=LocalEstimation_B(I, b_Label, Obj_Label, R, MaxIter);

b_Sigma=Cal_b_Sigma(b, I, b_Label, Obj_Label, R);

end
%================================================================================================


%================================================================================================
function [b b_Sigma]=LocalEstimation_A(I, b_Label, Obj_Label, R, MaxIter)
%
% O--------------------x
% |
% | 1: (x-R, y-R)           3: (x+R, y-R)
% |
% |                 (x, y)
% |
% | 2: (x-R, y+R)           4: (x+R, y+R)
% Y

[Ly Lx]=size(I);

[X Y]=meshgrid(-R:R, -R:R);
dist=X.^2+Y.^2;

Area_Thresh=9;

b=NaN(Ly, Lx);  
b_Sigma=NaN(Ly, Lx);

I_pad=padarray(I, [R R], 'replicate');

Im=I;
Im(b_Label==0)=0;
Im_pad=padarray(Im, [R R], 'replicate');

b_Label_pad=padarray(b_Label, [R R], 'replicate');

[YList XList]=find(Obj_Label);
TotalNum=length(XList);

y=YList(1); 
x=XList(1);
sub_I=I_pad(y:y+2*R, x:x+2*R);       
MeanValue=median(sub_I(:));
SigmaValue=sqrt(2*MeanValue); % theoritical value

for k=1:TotalNum
    y=YList(k); 
    x=XList(k);
    %----------------------------------------------------------------------------------
    % SegLabeled local mean estimation      
    % find the optimal W and b(y, x)
    sub_I=I_pad(y:y+2*R, x:x+2*R);           
    sub_b_Label=b_Label_pad(y:y+2*R, x:x+2*R);        
    subIdx=find(sub_b_Label>0); subIdx=subIdx(:);
    subNum=length(subIdx);        
    
    if subNum >= Area_Thresh
        Data=sub_I(subIdx);
        Weight=exp(-0.5*dist(subIdx)/(2*R)^2);        
        [MeanValue SigmaValue]=RobustLocalMean(Data, Weight, MaxIter);
                
    else
        sub_Im=Im_pad(y:y+2*R, x:x+2*R);
        subIdx=find(sub_Im>0); subIdx=subIdx(:);
        subNum=length(subIdx);   
        if subNum >= Area_Thresh
            Data=sub_Im(subIdx);
            Weight=exp(-0.5*dist(subIdx)/(2*R)^2);        
            [MeanValue SigmaValue]=RobustLocalMean(Data, Weight, MaxIter);
        end
    end
    
    Im(y+R, x+R)=MeanValue;
    
    b(y,x)=MeanValue;   
    b_Sigma(y,x)=SigmaValue;
    %----------------------------------------------------------------------------------                                    
end

end
%================================================================================================

%================================================================================================
function [b b_Sigma]=LocalEstimation_B(I, b_Label, Obj_Label, R, MaxIter)
%
% O--------------------x
% |
% | 1: (x-R, y-R)           3: (x+R, y-R)
% |
% |                 (x, y)
% |
% | 2: (x-R, y+R)           4: (x+R, y+R)
% Y
[Ly Lx]=size(I);
b=zeros(Ly, Lx);  
b_Sigma=NaN(Ly, Lx);

Im=I;
Im(b_Label==0)=0;

Area_Thresh=9;

[YList XList]=find(Obj_Label);
TotalNum=length(XList);
   
MeanValue=median(I(Obj_Label>0));
SigmaValue=sqrt(2*MeanValue); % theoritical value

for k=1:TotalNum
    y=YList(k); 
    x=XList(k);                                                                 
        
    xs=max(x-R, 1); xe=min(x+R, Lx);
    ys=max(y-R, 1); ye=min(y+R, Ly);   
              
    [sub_X sub_Y]=meshgrid(xs:xe, ys:ye);
        
    dist=(sub_X-x).^2+(sub_X-y).^2;                          
        
    %----------------------------------------------------------------------------------
    % find the optimal weight W  and bm2_t(y, x)                                  
    % weighted local linear regression                
    sub_I=I(ys:ye, xs:xe);
    sub_b_Label=b_Label(ys:ye, xs:xe);
    subIdx=find(sub_b_Label>0); subIdx=subIdx(:);
    subNum=length(subIdx);
    
    if subNum >= Area_Thresh
        Data=sub_I(subIdx);
        Weight=exp(-0.5*dist(subIdx)/(2*R)^2);
        [P SigmaValue]=RobustLocalPlaneFit(Data, sub_X(subIdx), sub_Y(subIdx), Weight, MaxIter);
        MeanValue=P(1)+P(2)*x+P(3)*y;
    else
        sub_Im=Im(ys:ye, xs:xe);        
        subIdx=find(sub_Im>0); subIdx=subIdx(:);
        subNum=length(subIdx); 
        if subNum >= Area_Thresh
            Data=sub_Im(subIdx);
            Weight=exp(-0.5*dist(subIdx)/(2*R)^2);
            [P SigmaValue]=RobustLocalPlaneFit(Data, sub_X(subIdx), sub_Y(subIdx), Weight, MaxIter);
            MeanValue=P(1)+P(2)*x+P(3)*y;
        end
        
    end
    
    Im(y,x)=MeanValue;
    
    b(y,x)=MeanValue;   
    b_Sigma(y,x)=SigmaValue;
    %----------------------------------------------------------------------------------                                                      

end

end
%================================================================================================


%================================================================================================
function [b b_Sigma]=LocalEstimation_C(I, b_Label, Obj_Label, R, MaxIter)
%
% O--------------------x
% |
% | 1: (x-R, y-R)           3: (x+R, y-R)
% |
% |                 (x, y)
% |
% | 2: (x-R, y+R)           4: (x+R, y+R)
% Y

[Ly Lx]=size(I);

[X Y]=meshgrid(-R:R, -R:R);
dist=X.^2+Y.^2;

b=NaN(Ly, Lx);  
b_Sigma=NaN(Ly, Lx);

I_pad=padarray(I, [R R], 'replicate');

b_Label_pad=padarray(b_Label, [R R], 'replicate');

[YList XList]=find(Obj_Label);
TotalNum=length(XList);

y=YList(1); 
x=XList(1);
sub_I=I_pad(y:y+2*R, x:x+2*R);       
MeanValue=median(sub_I(:));
SigmaValue=sqrt(2*MeanValue); % theoritical value

for k=1:TotalNum
    y=YList(k); 
    x=XList(k);
    %----------------------------------------------------------------------------------
    % SegLabeled local mean estimation      
    % find the optimal W and b(y, x)
    sub_I=I_pad(y:y+2*R, x:x+2*R);       
    sub_b_Label=b_Label_pad(y:y+2*R, x:x+2*R);        
    subIdx=find(sub_b_Label>0); subIdx=subIdx(:);
    subNum=length(subIdx);        
    
    if subNum > 0
        MeanValue=median(sub_I(subIdx));                
        SigmaValue=1.4826*median(abs(MeanValue-sub_I(subIdx)));
    end
    
    b(y,x)=MeanValue;   
    b_Sigma(y,x)=SigmaValue;
    %----------------------------------------------------------------------------------                                    
end

%------------------------------------------------------------------------------------------
%b_Sigma=Cal_b_Sigma(b, I, b_Label, Obj_Label, R);
%------------------------------------------------------------------------------------------
end
%================================================================================================


%================================================================================================
function b_Sigma=Cal_b_Sigma(b, I, b_Label, Obj_Label, R)
% calculate the sigma of b at (y, x), i.e., Median absolute deviation
[Ly Lx]=size(b);
b_Sigma=NaN(Ly, Lx);

[YList XList]=find(Obj_Label);
TotalNum=length(XList);

y=YList(1); 
x=XList(1);
SigmaValue=sqrt(2*b(y,x)); % theoritical value
    
for k=1:TotalNum
    y=YList(k); 
    x=XList(k);
                
    xs=max(x-R, 1); xe=min(x+R, Lx);
    ys=max(y-R, 1); ye=min(y+R, Ly);   
                
    sub_b=b(ys:ye, xs:xe);
    sub_I=I(ys:ye, xs:xe);      
    sub_b_Label=b_Label(ys:ye, xs:xe);
    subIdx=find(sub_b_Label>0);        
    subNum=length(subIdx);
            
    if subNum>0
        b_Sigma(y,x)=1.4826*median(abs(sub_b(subIdx)-sub_I(subIdx)));        
        SigmaValue=b_Sigma(y,x);
    else
        b_Sigma(y,x)=SigmaValue;
    end
end

end
%==========================================================================================================


%==========================================================================================================
function [MeanValue SigmaValue]=RobustLocalMean(Data, Weight, MaxIter)
 
W=ones(size(Data));    
W=W.*Weight+eps;
for Iteration=1:MaxIter

    W=W./sum(W);
    MeanValue=sum(W.*Data);
    SigmaValue=1.4826*median(abs(MeanValue-Data))+eps;            
    W=exp(-0.5*(Data-MeanValue).^2/SigmaValue^2).*Weight;        
    W=W+eps;    
end
end
%==========================================================================================================


%==========================================================================================================
function [P SigmaValue]=RobustLocalPlaneFit(Data, X, Y, Weight, MaxIter)

W=Weight;
Num=length(Data);
A=[ones(Num, 1) X Y];
for Iteration=1:MaxIter               
    P=lscov(A, Data, W);                              
    FitValue=P(1)+P(2)*X+P(3)*Y;  
    FitValue=max(FitValue, 0);          
    SigmaValue=1.4826*median(abs(Data-FitValue))+eps;
    W=exp(-0.5*(Data-FitValue).^2/SigmaValue^2).*Weight;
    W=W+eps;
end        
end
%==========================================================================================================