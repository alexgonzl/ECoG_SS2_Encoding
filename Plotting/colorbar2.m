function han = colorbar2(ax,CM)
% function creates a color bar with the colormap 2
% ax = axes
% CM = colormap

vert = 1;
Nlevels = size(CM,1);
axes(ax);

if vert
    ylim([1 Nlevels+1])
    xlim([0 1])
end
xLims = xlim;
yLims = ylim;

S = zeros(4,1);
if vert
    S(1) = xLims(1);
    S(2) = xLims(2);
    S(3) = xLims(2);
    S(4) = xLims(1);    
    X = S;
    Z = yLims;
else
    S(1) = yLims(1);
    S(2) = yLims(2);
    S(3) = yLims(2);
    S(4) = yLims(1);       
    Y = S;    
    Z = xLims;
end

cnt = 0;
for ii = 1:Nlevels
    if vert 
        Y(1) = Z(1)+cnt;
        Y(2) = Z(1)+cnt;
        Y(3) = Z(1)+cnt+1;
        Y(4) = Z(1)+cnt+1;
    else
        X(1) = Z(1)+cnt;
        X(2) = Z(1)+cnt;
        X(3) = Z(1)+cnt+1;
        X(4) = Z(1)+cnt+1;
    end
    
    p=patch(X,Y,CM(ii,:));
    p.EdgeColor='none';
    cnt = cnt + 1;
    
end
axis off

