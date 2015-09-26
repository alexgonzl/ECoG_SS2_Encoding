% Illustration of SPL attention
addpath Analysis/
addpath Plotting/
addpath lib/

dirPath         = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
fileName        = 'allERSPshgamGroupstimLocksublogPowernonLPCCh';

dataPath = [dirPath 'group/Spectral_Data/'];
load([dataPath fileName '.mat'])

t = data.trialTime;
xtick      = [0 0.4 0.8 1.2];
xticklabel = {'stim','0.4','0.8','1.2'};

%%
mu1 = 0.2;
sd1 = 0.12;
y1= gampdf(t,3,0.1);
y1=y1/max(y1);

mu1 = 0.35;
sd1 = 0.1;
y2 = gampdf(t,2.5,0.3);
y2 = y2/max(y2);

fP = '~/Google Drive/Research/ECoG_SS2e/';
figure(2);
set(gcf,'position',[100 200 750 100],'paperpositionmode','auto')
imagesc(t,1,y1);colormap hot
set(gca,'ytick',[],'xtick',info.xtick,'xticklabel',xticklabel)
set(gca,'fontsize',30,'FontWeight','bold')
print(gcf,'-dtiff',[fP 'SPL-perception'])

figure(3);
set(gcf,'position',[100 200 750 100],'paperpositionmode','auto')
imagesc(t,1,y2);colormap hot
set(gca,'ytick',[],'xtick',xtick,'xticklabel',xticklabel)
set(gca,'fontsize',30,'FontWeight','bold')
print(gcf,'-dtiff',[fP 'SPL-semantic'])


