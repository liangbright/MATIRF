function Weight=Module2_Cal_Weight(ParticleState, A_Data, b_Array)
    
[ParticleNum, AngleNum]=size(A_Data);
    
Weight=ones(ParticleNum, AngleNum);
    
for m=1:AngleNum
    
    Ib=b_Array(:,:,m);
    
    x=ParticleState(2,:);
    y=ParticleState(3,:);
    for k=1:ParticleNum
        Axy=A_Data(k,m);
        bxy=Ib(y(k),x(k));
            
        Weight(k,m)=Axy/(Axy+bxy);
    end    
end