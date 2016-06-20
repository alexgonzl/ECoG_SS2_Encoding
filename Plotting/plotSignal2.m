%%
function han=plotSignal2(t,x,cols)

nSigs = size(x,2); % each signal is a column

han=figure(); clf;
set(gcf,'paperpositionmode','auto','color','w')
ar = [600 150];
set(gcf,'paperUnits','points','papersize',ar,'paperposition',[0 0 ar])
set(gcf,'position',[200,300,ar]); hold on;
for ii = 1: nSigs
    plot(t,x(:,ii),'color',cols(ii,:),'linewidth',2)
end
set(gca,'linewidth',2,'fontsize',20)
xtick = round([t(1) mean(t) t(end)]*10)/10;
set(gca,'xtick',xtick)
axis tight
dx=(t(end)-t(1))*0.08;
xlim([t(1)-dx t(end)+dx])

yL = ylim;
dy=(abs(yL(2)-yL(1)))*0.08;
ylim([yL(1)-dy yL(end)+dy])

xlabel(' Time Epoch (s) ')
ylabel(' \mu V ')

end
