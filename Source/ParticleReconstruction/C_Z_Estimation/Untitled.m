for k=1:100
Sigma=0.1*k;    
x=-100:1:100;
y=-100:1:100;

h=exp(-(x.^2+y.^2)/Sigma^2);

sumh(k)=sum(h(:));
end

figure; plot(sumh); grid on