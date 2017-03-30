function Flag=IsNearBoundary(CellLabel, x, y, Radius)

[Ly Lx]=size(CellLabel);

% x1 ----x2
%  |  x  |
% x3 --- x4

x1=round(x-Radius); x2=round(x+Radius); 
y1=round(y-Radius); y3=round(y+Radius);

if x1<=0 || x2>Lx || y1<=0 || y3>Ly
    Flag=true;
else
    
    Label=CellLabel(y1:y3, x1:x2);
    
    if sum(Label==0) > 0
        Flag=true;
    else
        Flag=false;
    end
end