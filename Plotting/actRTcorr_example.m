

subjs = 1:4;

info = [];
info.cond1 	= cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
info.cond1 	= cellfun(@and,info.cond1,data.conds(subjs,5),'uniformoutput',0); % correctly remembered
info.cond2 	= cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
info.cond2 	= cellfun(@and,info.cond2,data.conds(subjs,5),'uniformoutput',0); % correct remembered
info.cond3  = cellfun(@or, info.cond1,info.cond2,'uniformoutput',0); % collapse across semantic task

info.xtick      = [0 0.4 0.8 1.2];
info.xticklabel = {'stim','0.4','0.8','1.2'};
info.Bins       = data.Bins;
info.toi        = [0 1.5];
info.groupNames = {'IPS','SPL','AG'};
info.cols       = 'k';

info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';

%% plot a single channel

s  = 3;
ch = 12;
tr = info.cond3{s};
prefix = sprintf('s%ich%i_%s%s_',s,ch,'stim','hgam');

info.fileName   = sprintf('%s_singleChanTC',prefix);
info.yLimits    = [-1 5];
X = [];
X{1} = squeeze(data.ERP{s}(ch,tr,:));
t = data.trialTime;
info.xlabel= ' Time(s) ';
info.ylabel= ' HGP(dB) ';
plotWrapper(X,t,info)


%%
rts=data.studyRTs{s}(tr);
[h,i]=hist((rts),20);

figure(2);
set(gcf,'paperpositionMode','auto')
h =bar(i,h);
set(h,'barwidth',1, 'edgecolor','none','faceColor','k')
set(gca,'fontSize',16,'ytick',[0:4:12],'xtick',0:1:5)
xlabel(' RTs(s) ')
ylabel(' Trial Counts ')

print(gcf, '-depsc2', [info.savePath prefix 'studyRTshist'])

%%

X = data.BinERP{s}(ch,tr,7)';
figure(3);
h = scatter(rts,X);
set(h,'markerfaceColor','k','markeredgeColor','none','sizeData',40)
set(gca,'fontSize',18,'ytick',-6:3:9,'xtick',0:1:5)
xlabel(' RTs(s) ')
ylabel(' HGP(dB) ')
l=lsline;
set(l,'linewidth',3)
r = corr(log(rts),X);
text(3.5,7,sprintf('r = %0.2f',r),'fontsize',16)

print(gcf, '-depsc2', [info.savePath prefix 'scatterCorrPlot'])
%%

t = mean(data.Bins,2);
X = corr(squeeze(data.BinERP{s}(ch,tr,:)),rts);

figure(4); clf;hold on;
plot(data.trialDur,[0 0],'--', 'linewidth',2,'color',0.3*ones(1,3))
plot([0 0],[-.1 .1],'--', 'linewidth',2,'color',0.3*ones(1,3))
h=plot(t,X);
set(h,'lineWidth',4,'color','k')
set(gca,'fontSize',16,'ytick',[-0.4 -0.2 0 0.2 0.4],'xtick',info.xtick,'xticklabel',info.xticklabel)
xlim(data.trialDur)
xlabel( ' Time (s) ')
ylabel( ' Corr (r) ')

print(gcf, '-depsc2', [info.savePath prefix 'CorrTC'])



