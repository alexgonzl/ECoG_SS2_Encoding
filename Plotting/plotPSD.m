 
function han = plotPSD(f,pxx)
han=figure(); clf;
set(gcf,'paperpositionmode','auto','color','w')
ar = [400 300];
set(gcf,'paperUnits','points','papersize',ar,'paperposition',[0 0 ar])
set(gcf,'position',[200,500,ar]); hold on;
plot(f,10*log10(pxx),'k','linewidth',2)
set(gca,'linewidth',2,'fontsize',20)
xtick = round(linspace(f(1),f(end),5));
set(gca,'xtick',xtick,'ytick',-80:20:80)
xlabel(' Frequency (Hz) ')
ylabel(' Power (dB) ')
axis tight
end
