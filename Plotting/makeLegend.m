function han=makeLegend(colors,textStr)

figure(); clf;

set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[200 200],'paperposition',[0 0 200 200])
set(gcf,'position',[200, 200, 200, 200])
set(gca,'units','points','position',[50 50 100 100]); hold on;

nLabels = size(colors,1);


cnt = 1;
for ii = nLabels:-1:1
    s = scatter(1,cnt,'o');
    s.MarkerFaceAlpha   = 0.9;
    s.MarkerEdgeAlpha   = 0.9;
    s.SizeData          = 120;
    s.MarkerEdgeColor   = colors(ii,:);
    s.MarkerFaceColor   = colors(ii,:);
    
    t=text(3,cnt,textStr{ii},'fontsize',20);
    cnt = cnt+1;
end
xlim([0 10])
ylim([0,nLabels+1])
axis off
han = gca;