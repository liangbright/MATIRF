Depth=[150 200 300 400 500];
for k=1:5
    d=Depth(k);
    R=50;
    y1(k)=((R/d)*cosh(R/d)-sinh(R/d))*(d/R)^3;
    R=60;
    y2(k)=((R/d)*cosh(R/d)-sinh(R/d))*(d/R)^3;
    
    R=100;
    y3(k)=((R/d)*cosh(R/d)-sinh(R/d))*(d/R)^3;
end

m=1:k;
figure; plot(m, y1, 'r', m, y2, 'b', m, y3, 'g'); grid on