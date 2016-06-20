function GLM_PCA_SB_CompPlots(opts)

% load data and PCA results
fileName = ['allMBAnalysis' opts.lock 'sublogPowernonLPCch'];
load([opts.dataPath opts.lock '/' fileName '.mat'])
fileName = ['PCATrialDecomp-SBAnalysis' opts.lock 'sublogPowernonLPCch'];
load([opts.dataPath opts.lock '/' fileName '.mat'])
pcadat = out; clear out;
% load electrode locations
load('~/Google Drive/Research/ECoG_SS2e/data_results/Renderings/electrodeLocs.mat')

%rThr    = opts.rThr;
freqs   = 1:data.nBands;
time    = mean(data.Bins(data.AnalysisBins,:),2);
nFeat   = pcadat.nFeat;
nChans  = pcadat.nChans;
nComps  = pcadat.nComps;
nBins   = data.nBins;
nBands  = data.nBands;
rois    = data.ROIid;
nROIs   = numel(unique(rois));
% colormap settings
clev    = [-4:-1 1:4];
y       =     brewermap(numel(clev)-1,'*RdBu');
rcolmap = brewermap(1000,'*RdBu');
nBins   = sum(data.AnalysisBins);
ROIcolors = [240 35 17; 2 93 140;122 179 23]/255;

%% bar plot by region
if 1
opts2               = [];
opts2.colors        = ROIcolors;
opts2.aspectRatio   = [500 300];
opts2.xTick         = 1:6;
opts2.xticksLabels  = {'\delta','\theta','\alpha','\beta','l-\gamma','h-\gamma'};
%opts2.xTickLabel    = {'IPS','SPL','AG'};
opts2.xLimits       = [0.4 6.6];
opts2.yLabel        = ' R^2 ';
opts2.yLimits       = [0 20];
opts2.yTicks        = [ 5 10 15];
main_xlocks         = [0.75,1,1.25];

p4aName = ['StudyGLM_PCA_SB_R2_ROIs'];
p4bName = ['TestGLM_PCA_SB_R2_ROIs'];
nROIs = 3;
x =  pcadat.StudyGLMsChanR2*100;
M = cell(nROIs,nBands);

figure(1); clf;hold on;
set(gcf,'paperpositionmode','auto','position',[100 100 opts2.aspectRatio])
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[opts2.aspectRatio],'paperposition',[0 0 opts2.aspectRatio])
set(gcf,'position',[100,100,opts2.aspectRatio]); 

for rr =1:nROIs
    for bb = 1:nBands
        M{rr,bb} = x(rois==rr,bb);
        N{rr,bb} = sum(~isnan(M{rr,bb}));
        [me(rr,bb),se(rr,bb)] = grpstats(M{rr,bb},ones(numel(M{rr,bb}),1),{'mean','sem'});
        
        bar(main_xlocks(rr)+bb-1,me(rr,bb),0.25,'FaceColor',opts2.colors(rr,:),'edgeColor', 'none', 'basevalue',0,'ShowBaseLine','off')
        plot([1 1]*main_xlocks(rr)+bb-1, [-se(rr,bb) se(rr,bb)]+me(rr,bb), 'color',[0 0 0],'linewidth',4)
    end
end
set(gca,'xTick',opts2.xTick,'xTickLabel',opts2.xticksLabels,'YTick',opts2.yTicks,'fontsize',24)
set(gca,'linewidth',2)
axis tight
xlabel('Band')
ylabel('R^2')
print(gcf,'-dpdf', [opts.savePath p4aName])

x =  pcadat.TestGLMsChanR2*100;
figure(1); clf;hold on;
set(gcf,'paperpositionmode','auto','position',[100 100 opts2.aspectRatio])
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[opts2.aspectRatio],'paperposition',[0 0 opts2.aspectRatio])
set(gcf,'position',[100,100,opts2.aspectRatio]); 

for rr = 1:nROIs
    for bb = 1:nBands
        M{rr,bb} = x(rois==rr,bb);
        N{rr,bb} = sum(~isnan(M{rr,bb}));
        [me(rr,bb),se(rr,bb)] = grpstats(M{rr,bb},ones(numel(M{rr,bb}),1),{'mean','sem'});
        
        bar(main_xlocks(rr)+bb-1,me(rr,bb),0.25,'FaceColor',opts2.colors(rr,:),'edgeColor', 'none', 'basevalue',0,'ShowBaseLine','off')
        plot([1 1]*main_xlocks(rr)+bb-1, [-se(rr,bb) se(rr,bb)]+me(rr,bb), 'color',[0 0 0],'linewidth',4)
    end
end
set(gca,'xTick',opts2.xTick,'xTickLabel',opts2.xticksLabels,'YTick',opts2.yTicks,'fontsize',24)
set(gca,'linewidth',2)
axis tight
xlabel('Band')
ylabel('R^2')
print(gcf,'-dpdf', [opts.savePath p4bName])
end