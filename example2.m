load('/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/data_results/preStim2/allMBAnalysispreStim2sublogPowernonLPCch.mat')
load('/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/data_results/preStim2/PCATrialDecomp-MBAnalysis2_KmeanspreStim2sublogPowernonLPCch.mat')

%load('/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/data_results/stim/allMBAnalysisstimsublogPowernonLPCch.mat')
%load('/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/data_results/stim/PCATrialDecomp-MBAnalysis2_KmeansstimsublogPowernonLPCch.mat')

%% find a large correlation

[ch,co]=find(out.CorrTestRTs>0.4);
cc =5;
ch = ch(cc);co=co(cc);

ss = data.subjChans(ch);
subjChan = find(find(data.subjChans==ss)==ch);
rts = out.tRTs{ss};
Comp = squeeze(out.Comps{ss}(subjChan,:,co))';

%% 
ar = [300 300];
han=figure(1); clf;
set(gcf,'paperpositionmode','auto','color','w')
set(gcf,'paperUnits','points','papersize',ar,'paperposition',[0 0 ar])
set(gcf,'position',[200,300,ar]); hold on;
s = scatter(Comp,-log10(rts));
s.MarkerFaceColor = [0.4 0.4 0.4];
s.MarkerEdgeColor = 'none';
s.SizeData = 150;

axis square

set(gca,'linewidth',2,'fontsize',20)
xlabel(' Prin Comp ')
ylabel(' -log_{10} (RT) ')
axis tight


xLim = xlim;
yLim = ylim;
dx = xLim(2)-xLim(1);
dy = yLim(2)-yLim(1);

xlim([xLim(1)-dx*0.1 xLim(2)+dx*0.1])
ylim([yLim(1)-dy*0.1 yLim(2)+dy*0.1])

fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(han,'-dpdf', [fPath 'Scatter_Comp_TestRT'])

%%
freqs=1:6;
figure();
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[500 300],'paperposition',[0 0 500 300])
set(gcf,'position',[100,100,500,300]); % 100pt all around margin
a   = axes('position',[0.15 0.15 0.70 0.5]);

X=reshape(squeeze(out.Projections(ch,:,co)),[6,out.nFeat/6]);
t=mean(data.Bins,2); t=t(data.AnalysisBins);
xtick      = [-.8 -0.4 0 0.4 0.8 1.2];
xticklabel = {'-0.8','-0.4','stim','0.4','0.8','1.2'};
%xtick      = [0 0.4 0.8 1.2];
%xticklabel = {'stim','0.4','0.8','1.2'};

clev    = [-4:-1 1:4];
clev2 = quantile(X(:),(clev+5)/10);
cmm       =     brewermap(numel(clev)-1,'*RdBu');
rcolmap = brewermap(1000,'*RdBu');
h=contourfcmap(t,freqs, X,clev2,cmm,'lo',rcolmap(1,:),'hi', rcolmap(end,:),'method','calccontour');
for i = 1:numel(h.l)
    h.l(i).Color='none';
end
hold on;
%plot([0 0],[1 6],'k','linewidth',4)
yticksLabels   = {'\delta','\theta','\alpha','\beta','l-\gamma','h-\gamma'};
set(gca,'fontsize',20,'ytick',freqs,'ytickLabel',yticksLabels,'xtick',xtick,'xticklabel',xticklabel)
plot([0 0],[1 6],'k','linewidth',4)
print(gcf,'-dtiff','-r300', [fPath 'Spectrogram_SelComp_TestRT2'])

%%

t = linspace(0,1,2000);
x1 = cos(2*pi*10*t);
x2 = cos(2*pi*6*t);

figure();
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 150],'paperposition',[0 0 600 150])
set(gcf,'position',[100,100,600,150]); 

plot(t,x1,'linewidth',10,'color', 'k')
set(gca,'xtick',[0:0.5:2],'fontsize',20,'linewidth',2,'ycolor','w','box','off')
xlabel( 'Time (s) ')
axis off

print(gcf,'-dpdf', [fPath 'AlphaRhythm'])

cla;
plot(t,x2,'linewidth',10,'color','k')
set(gca,'xtick',[0:0.5:2],'fontsize',20,'linewidth',2,'ycolor','w','box','off')
xlabel( 'Time (s) ')
axis off

print(gcf,'-dpdf', [fPath 'ThetaRhythm'])



