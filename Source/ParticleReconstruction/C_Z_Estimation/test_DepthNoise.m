Depth=[100 200 300 400 500];
d=Depth+sqrt(Depth).*randn(size(Depth));

%%
Depth=100;
counter=0;
std_max=[];
for x=0:0.001:0.002
for z=1:1000
    d=Depth+(Depth*x).*randn(10000,1);        
    y=exp(-z./d);
    y_std(z)=std(y);
end     
counter=counter+1;
std_max(counter)=max(y_std);
end
%%
figure; plot(y_std); grid on
figure; plot([0:0.001:0.002], std_max); grid on
%%