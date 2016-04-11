function figHandle = plotWrapper(X,t,info)

inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';


% defaults
if ~isfield(info,'centerType')
    info.centerType 	= 'mean';
end
if ~isfield(info,'yLimits')
    info.yLimits        = [];
end
if ~isfield(info,'yCenter')
    info.yCenter        = 0;
end
if ~isfield(info,'yRefLimits')
    info.yRefLimits     = [];
end
if ~isfield(info, 'smoother')
    info.smoother		= 'loess';
    info.smootherSpan 	= 0.15;
end
if ~isfield(info,'errorBars')
    info.errorBars      = 'se';
end
if ~isfield(info,'cols')
    info.cols           = 'rbgk';
end
if ~isfield(info,'savePath')
    info.savePath       = pwd;
elseif ~exist(info.savePath,'dir')
    mkdir(info.savePath);
end


figHandle = figure(1); clf;
set(gcf,'position',[100 200 750 500],'paperpositionmode','auto')
h = plotNTraces(X,t,'cols',info.cols,...
    'smoother',info.smoother,'smootherSpan',info.smootherSpan,...
    'centerType',info.centerType,'errorBars',info.errorBars,...
    'yLimits',info.yLimits);

mainLines = zeros(numel(X),1);
for ii = 1:numel(X)
    mainLines(ii) = h.(['h' num2str(ii)]).mainLine;
end

if isfield(info,'PVals') && isfield(info,'Bins') && isfield(info,'alpha')
    p 	= info.PVals;
    sigBins = info.Bins(p<info.alpha,:);
    if isfield(info,'toi')
        sigBins = sigBins([sigBins(:,1)>=info.toi(1)&sigBins(:,2)<=info.toi(2)],:);
    end
    hold on;
    ylims = ylim;
    for ii = 1:size(sigBins,1)
        plot(sigBins(ii,:),ones(1,2)*(ylims(2)-ylims(2)*0.1),'linewidth',2,'color',0.3*ones(3,1))
        plot(mean(sigBins(ii,:)),(ylims(2)-ylims(2)*0.1),'*','linewidth',2,'color',0.1*ones(3,1))
        
    end
end

if isfield(info,'yAxisRightLoc')
    if info.yAxisRightLoc
        set(gca,'YAxisLocation','right')
    end
end

if isfield(info,'xtick') && isfield(info,'xticklabel')
	set(gca,'xtick',info.xtick)
	set(gca,'xticklabel',info.xticklabel)
end

if isfield(info,'legend')
    l = legend(mainLines, info.legend,'location','best');
    set(l,'box','off')
end

if isfield(info,'xlabel')
    axes('position',[0.4 0 0.2 0.1])
    text(0.2,0.3, info.xlabel,'fontsize',22)
    axis off
end
if isfield(info,'ylabel')
    axes('position',[0 0.4 0.1 0.2])
    text(0.45,0.1, info.ylabel,'fontsize',22,'rotation',90)
    axis off
   % ylabel(info.ylabel)
end

set(gca,'fontSize',24)

cPath = pwd;
cd(info.savePath)
addpath(cPath)
addpath([cPath '/Plotting/'])

try
    fN = strcat(info.fileName,'.svg');
    plot2svg(fN,h.f)
    
    eval(['!' inkscapePath ' -z ' fN ' --export-pdf=' info.fileName '.pdf'])
    eval(['!rm ' fN])
catch
end
cd(cPath)