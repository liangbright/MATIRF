d1=100;
d2=400;
value=[];
C_g=0;
z_g=0;
value_min=inf;
for C=100
    for z=0:10:1000
        value(end+1)=sqrt((0-C*exp(-z/d1))^2+(100-C*exp(-z/d2))^2);
        
        if value_min>value(end)
            value_min=value(end);
            C_g=C;
            z_g=z;
        end
    end
end
%%
figure; plot(value,'o-');