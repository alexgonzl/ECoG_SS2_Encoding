function behavPlots_prelim(savePath)
%% Plot the results by subject.

dataPath  = '~/Google Drive/Research/ECoG_SS2e/Behavior/summary/';
load([dataPath 'behavSummary.mat'])

inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';
%fP = '~/Google Drive/Research/ECoG_SS2e/Plots/Behavior/';
fP = savePath;

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

% Figure: combined plots

% axis 1, accuracy.
y = AC;
opts = [];
opts.shapes = shapes;
opts.f_han = figure(1); clf;
set(gcf, 'position', [100 100 800 400],'paperpositionmode','auto');
opts.markerSize = 125;
opts.ax_han =  axes('position',[0.05 0.1 0.15 0.85]);
opts.yLims = [-0.1 1.1];
opts.yticks = [0 0.5 1];
opts.text = sprintf(' M = %0.2f',mean(y));
%opts.text_loc = [0.5,0.2];
opts.xLabel = ' ACC ';
opts.yLabel = '';
opts.xticks = '';
opts.inkspace = inkscapePath;

coreplot(y,opts);

%---%
% axis 2, RTs (encoding)
y = RTe;
opts.ax_han = axes('position',[0.30 0.1 0.15 0.85]); hold on;
opts.yLims = [-0.2 2.6];
opts.xLabel = ' RTs-enc (s) ';
opts.yticks = [0 1 2];
opts.text = sprintf(' M = %0.2fs',mean(RTe));
%opts.text_loc = [0.5 0.2];

coreplot(y,opts);

% RTs (retrieval)
y = RTr;
opts.ax_han = axes('position',[0.55 0.1 0.15 0.85]); hold on;
opts.xLabel = ' RTs-ret (s) ';
opts.text = sprintf(' M = %0.2fs',mean(RTr));

coreplot(y,opts);

% RT-corr ( enc-ret) and save
y = RTsCorr;
opts.ax_han = axes('position',[0.80 0.1 0.15 0.85]); hold on;
opts.yLims = [-0.6 0.6];
opts.xLabel = ' enc-ret (r) ';
opts.text = sprintf(' M = %0.2f',mean(RTsCorr));
%opts.text_loc = [0.5,-.4];
opts.yticks = [-0.5 0 0.5];


opts.filePath = fP;
opts.fN =  'encodingPerfScatter';
coreplot(y,opts);

%% individual plots
% accuracy
y = AC;
opts = [];
opts.shapes = shapes;
opts.f_han = figure(1); clf;
set(gcf, 'position', [100 100 500 400],'paperpositionmode','auto');
opts.markerSize = 200;
opts.ax_han = axes('position',[0.1 0.1 0.38 0.85]); hold on;

opts.yLims = [-0.1 1.1];
opts.yticks = [0 0.5 1];
opts.text = sprintf(' M = %0.2f',mean(y));
opts.xLabel = ' ACC ';
opts.yLabel = '';
opts.xticks = '';
opts.inkspace = inkscapePath;

coreplot(y,opts);

% RTs (encoding)
y = RTe;
opts.ax_han = axes('position',[0.55 0.1 0.38 0.85]); hold on;
opts.yLims = [-0.2 2.2];
opts.yticks = [0 1 2];
opts.xLabel = ' RTs-enc (s) ';
opts.text = sprintf(' M = %0.2f',mean(y));
opts.filePath = fP;
opts.fN =  'encodingPerfScatter_1';

coreplot(y,opts);

%%
% RTs (retrieval)
y = RTr;
opts = [];
opts.shapes = shapes;
opts.markerSize = 200;
opts.f_han = figure(2); clf;
set(gcf, 'position', [100 100 500 400],'paperpositionmode','auto');
opts.ax_han = axes('position',[0.1 0.1 0.38 0.85]); hold on;
opts.yLims = [-0.2 2.6];
opts.yticks = [0 1 2];
opts.xLabel =' RTs-ret (s) ';
opts.yLabel = '';
opts.xticks = '';
opts.text = sprintf(' M = %0.2f',mean(y));

coreplot(y,opts);

% RT-corr ( enc-ret)
y = RTsCorr;
opts.ax_han = axes('position',[0.55 0.1 0.38 0.85]); hold on;
opts.yLims = [-0.6 0.6];
opts.yticks = [-0.5 0 0.5];
opts.xLabel =' enc-ret (r) ';
opts.text = sprintf(' M = %0.2f',mean(y));
opts.filePath = fP;
opts.fN =  'encodingPerfScatter_2';
opts.inkspace = inkscapePath;

coreplot(y,opts);

end

%%---------------------------------------------------------------------------%%
function [f_han, ax_han] = coreplot(y,opts)
% internal core plotting function for behavioral inputs
% this only meant for subject level data, where x position is irrelevant.

if isfield(opts,'markerSize')
    markerSize = opts.markerSize;
else
    markerSize = 200;
end

if isfield(opts,'shapes')
    subj_shapes = opts.shapes;
else
    subj_shapes = repmat('o',1,10);
end

if isfield(opts,'f_han')
    f_han = opts.f_han;
else
    f_han = figure(1); clf;
    set(f_han, 'position', [100 100 500 400],'paperpositionmode','auto');
end

if isfield(opts , 'ax_han')
    ax_han = opts.ax_han; hold on;
else
    ax_han = axes('position',[0.55 0.1 0.50 0.85]); hold on;
end

nSubjs = numel(y);
xpos = (1:nSubjs)/nSubjs*0.75;

% actual plotting
figure(f_han);
axes(ax_han);
xlim([0 1]);
plot([0.05 0.95],[1 1]*mean(y),'linewidth',3,'color',0.4*ones(1,3))
for ss= 1:nSubjs
    s=scatter(xpos(ss),y(ss),subj_shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',markerSize, 'lineWidth',3)
end

if isfield(opts,'yLims')
    ylim(opts.yLims);
    yLims = ylim;
else
end

if isfield(opts, 'xLabel')
    xlabel(opts.xLabel);
end
if isfield(opts, 'yLabel')
    ylabel(opts.yLabel);
end

if isfield(opts,'text')
    text_str = opts.text;
    if isfield(opts,'text_loc')
        text_loc = opts.text_loc;
    else
        text_loc = [0 0];
        text_loc(1) = 0.5; % middle of the panel
        text_loc(2) = yLims(1)+range(yLims)*0.2; % 20% from the bottom
    end

    text(text_loc(1),text_loc(2),text_str, 'fontsize',20,'horizontalalignment','center');
end

if isfield(opts,'yticks')
    set(ax_han,'ytick',opts.yticks);
end
   
if isfield(opts,'xticks')
    set(ax_han,'xtick',opts.xticks);
end

set(ax_han,'fontsize',20,'linewidth',2);

% save plot
if isfield(opts,'fN')
    fN = opts.fN;
    cPath = pwd;
    cd(opts.filePath)
    addpath(cPath)
    addpath([cPath '/Plotting/'])

    print(gcf,'-dsvg',fN)
    eval(['!' opts.inkspace ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
    eval(['!rm ' fN '.svg'])

    cd(cPath)
end

return % end coreplot
end
