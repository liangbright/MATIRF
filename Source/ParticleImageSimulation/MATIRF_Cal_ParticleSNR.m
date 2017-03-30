function [SNR, SNRPerAngle]=MATIRF_Cal_ParticleSNR(C, z, R, BgDensity, CellThickness, MATIRF_Param)

Depth=MATIRF_Param.Depth;
AngleNum=length(Depth);
I0=MATIRF_Param.I0;
SNRPerAngle=zeros(AngleNum,1);
for m=1:AngleNum               
        
    d=Depth(m);        
    
    if ~isempty(R)
        Alpha=3*((R/d).*cosh(R/d)-sinh(R/d)).*(d./R).^3;    
    else
        Alpha=1;
    end
    
    Am=I0(m)*C.*Alpha.*exp(-z/d);   
    
    bm=I0(m)*min(Depth(m),CellThickness)*BgDensity;
   
    SNRPerAngle(m)=mean(Am./sqrt(2*(Am+bm)));
end
SNR=mean(SNRPerAngle);