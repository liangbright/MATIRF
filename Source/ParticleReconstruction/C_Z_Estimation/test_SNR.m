bxy=40;
Axy=0:100;
        
y1=1-exp(-Axy.^2./(2*(bxy)));

y2=1-exp(-Axy.^2./(2*(Axy+bxy)));

y3=Axy./(Axy+bxy);


figure; plot(Axy, y1, 'g', Axy, y2, 'r', Axy, y3, 'b'); grid on