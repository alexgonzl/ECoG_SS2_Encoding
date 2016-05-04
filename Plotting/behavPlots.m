%% Plot the results by subject.
function behavPlots()

dataPath  = '~/Google Drive/Research/ECoG_SS2e/Behavior/summary/';
load([dataPath 'behavSummary.mat'])

inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';
fP = '~/Google Drive/Research/ECoG_SS2e/Plots/Behavior/';

nSubjs=numel(behavPerf.subjects);
AC = zeros(nSubjs,1);
RTe = zeros(nSubjs,1);
RTr = zeros(nSubjs,1);
RTsCorr = behavPerf.studytestRTsCorr(1:nSubjs);
for ss = 1:nSubjs
    AC(ss) = behavPerf.SubjStudyPerf(ss).AC;
    RTe(ss) = mean([behavPerf.Subj_studyRT_EncPerf(ss,:).mean]);
    RTr(ss) = mean([behavPerf.Subj_testRT_EncPerf(ss,:).mean]);
end

shapes = 'ods><x+^';
figure(1); clf; set(gcf, 'position', [100 100 800 400],'paperpositionmode','auto');
sD = 200;
% accuracy plot
axes('position',[0.05 0.1 0.18 0.85]); hold on;
xlim([0 1]); ylim([-0.1 1.1])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,AC(ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
end
plot([0.2 0.8],[1 1]*mean(AC),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 0.5 1],'linewidth',2);
xlabel(' ACC ')
text(0.5,0.1,sprintf(' M = %0.2f',mean(AC)),'fontsize',20,'horizontalalignment','center')


% RTs (encoding)
axes('position',[0.28 0.1 0.18 0.85]); hold on;
xlim([0 1]); ylim([-0.2 2.6])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,RTe(ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
end
plot([0.2 0.8],[1 1]*mean(RTe),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 1 2],'linewidth',2);
xlabel(' RTs-enc (s) ')
text(0.5,0.2,sprintf(' M = %0.2fs',mean(RTe)),'fontsize',20,'horizontalalignment','center')

% RTs (retrieval)
axes('position',[0.51 0.1 0.18 0.85]); hold on;
xlim([0 1]); ylim([-0.2 2.6])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,RTr(ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
end
plot([0.2 0.8],[1 1]*mean(RTr),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 1 2],'linewidth',2);
xlabel(' RTs-ret (s) ')
text(0.5,0.2,sprintf(' M = %0.2fs',mean(RTr)),'fontsize',20,'horizontalalignment','center')

% RT-corr ( enc-ret)
axes('position',[0.74 0.1 0.18 0.85]); hold on;
xlim([0 1]); ylim([-0.6 0.6])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,RTsCorr(ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
end
plot([0.2 0.8],[1 1]*mean(RTsCorr),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[-0.5 0 0.5],'linewidth',2);
xlabel(' enc-ret (r) ')
text(0.5,-.4,sprintf(' M = %0.2f',mean(RTsCorr)),'fontsize',20,'horizontalalignment','center')


fN = 'encodingPerfScatter';
cPath = pwd;
cd(fP)
addpath(cPath)
addpath([cPath '/Plotting/'])

print(gcf,'-dsvg',fN)
eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
eval(['!rm ' fN '.svg'])

cd(cPath)

%%

figure(1); clf; set(gcf, 'position', [100 100 500 400],'paperpositionmode','auto');
sD = 200;
% accuracy plot
axes('position',[0.1 0.1 0.40 0.85]); hold on;
xlim([0 1]); ylim([-0.1 1.1])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,AC(ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
end
plot([0.2 0.8],[1 1]*mean(AC),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 0.5 1],'linewidth',2);
xlabel(' ACC ')
text(0.5,0.1,sprintf(' M = %0.2f',mean(AC)),'fontsize',20,'horizontalalignment','center')


% RTs (encoding)
axes('position',[0.55 0.1 0.50 0.85]); hold on;
xlim([0 1]); ylim([-0.2 2.2])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,RTe(ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
end
plot([0.2 0.8],[1 1]*mean(RTe),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 1 2],'linewidth',2);
xlabel(' RTs-enc (s) ')
text(0.5,0.2,sprintf(' M = %0.2fs',mean(RTe)),'fontsize',20,'horizontalalignment','center')


fN = 'encodingPerfScatter_1';
cPath = pwd;
cd(fP)
addpath(cPath)
addpath([cPath '/Plotting/'])

print(gcf,'-dsvg',fN)
eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
eval(['!rm ' fN '.svg'])

cd(cPath)

figure(2); clf; set(gcf, 'position', [100 100 500 400],'paperpositionmode','auto');
sD = 200;
% RTs (retrieval)
axes('position',[0.1 0.1 0.40 0.85]); hold on;
xlim([0 1]); ylim([-0.2 2.6])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,RTr(ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
end
plot([0.2 0.8],[1 1]*mean(RTr),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 1 2],'linewidth',2);
xlabel(' RTs-ret (s) ')
text(0.5,0.2,sprintf(' M = %0.2fs',mean(RTr)),'fontsize',20,'horizontalalignment','center')

% RT-corr ( enc-ret)
axes('position',[0.55 0.1 0.50 0.85]); hold on;
xlim([0 1]); ylim([-0.6 0.6])
for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,RTsCorr(ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
end
plot([0.2 0.8],[1 1]*mean(RTsCorr),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[-0.5 0 0.5],'linewidth',2);
xlabel(' enc-ret (r) ')
text(0.5,-.4,sprintf(' M = %0.2f',mean(RTsCorr)),'fontsize',20,'horizontalalignment','center')

fN = 'encodingPerfScatter_2';
cPath = pwd;
cd(fP)
addpath(cPath)
addpath([cPath '/Plotting/'])

print(gcf,'-dsvg',fN)
eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
eval(['!rm ' fN '.svg'])

cd(cPath)