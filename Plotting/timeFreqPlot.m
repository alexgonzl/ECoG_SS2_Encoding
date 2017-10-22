function timeFreqPlot(time,freqs, X, clev,lock)

switch lock
    case 'stim'
        xtick      = [0 0.4 0.8 1.2];
        xticklabel = {'stim','0.4','0.8','1.2'};
    case 'RT'
        yAxisRightLoc = 1;
        xtick      = [-0.8 -0.4 0 0.4];
        xticklabel = {'-0.8','-0.4','resp','0.4'};
    case  'preStim2'
        xtick      = [-.8 -0.4 0 0.4 0.8 1.2];
        xticklabel = {'-0.8','-0.4','stim','0.4','0.8','1.2'};
    case  'preStim'
        xtick      = [-1.2 -.8 -0.4 0 0.4 0.8 1.2];        
        xticklabel = {'-1.2','-0.8','-0.4','stim','0.4','0.8','1.2'};
end

yticksLabels   = {'\delta','\theta','\alpha','\beta','l-\gamma','h-\gamma'};

y = brewermap(numel(clev)-1,'*RdBu');
y2 = brewermap(1000,'*RdBu');

axes('position',[0.1 0.1 0.85 0.8])
h=contourfcmap(time,freqs,X,clev,y, 'lo',y2(1,:),'hi', y2(end,:),'method','calccontour');
for i = 1:numel(h.l)
    h.l(i).Color='none';
end

hold on;
plot([0 0],[1 6],'k','linewidth',4)
set(gca,'fontsize',20,'YTick',freqs,'ytickLabel',yticksLabels)
set(gca,'XTick',xtick,'xtickLabel',xticklabel)

