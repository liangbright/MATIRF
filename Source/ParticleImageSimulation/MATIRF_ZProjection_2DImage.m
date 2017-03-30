function I=MATIRF_ZProjection_2DImage(Lx, Ly, IdxList_x, IdxList_y, xList, yList, z_min, z_max, ParticleFeature, MATIRF_Param, AngleIndex)

I=zeros(Ly,Lx);
for k=1:length(xList)        
    for n=1:length(yList)                                   
        Value=MATIRF_ZProjection_SinglePixel(xList(k), yList(n), ...
                                             z_min, z_max, ParticleFeature, MATIRF_Param, AngleIndex);    
        I(IdxList_y(n), IdxList_x(k))=Value;                                 
    end    
end