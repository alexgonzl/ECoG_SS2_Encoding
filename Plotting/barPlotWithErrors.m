function ha = barPlotWithErrors(M,opts)
% ha = barPlotWithErrors(M)
%
% M -> cell array with the data. Each cell entry becomes a column in the
% bargraph
% opts -> options:
%   aspectRatio
%   colors
%   baseLine
%   yLimits
%   xPositions

inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';

nBars   = numel(M);
m       = zeros(nBars,1);
se      = zeros(nBars,1);

if ~isfield(opts,'aspectRatio')
    opts.aspectRatio = [200 600];
end

if ~isfield(opts,'colors')
    opts.colors = zeros(nBars,3);
end

if ~isfield(opts,'baseLine')
    opts.baseLine =0;
end

if ~isfield(opts,'xPositions')
    opts.xPositions = 1:nBars;
end
ha=figure(); clf;hold on;
set(gcf,'paperpositionmode','auto','position',[100 100 opts.aspectRatio])

for ba = 1:nBars
    
    X = M{ba};
    n = size(X,1);
    [m(ba) se(ba)] = grpstats(X,ones(n,1),{'mean','sem'});
    
    bar(opts.xPositions(ba),m(ba),1,'FaceColor',opts.colors(ba,:),'edgeColor', 'none', 'basevalue',opts.baseLine,'ShowBaseLine','off')
    plot([1 1]*opts.xPositions(ba), [-se(ba) se(ba)]+m(ba), 'color',[0 0 0],'linewidth',4)
    
end

if isfield(opts,'yLimits')
    ylim(opts.yLimits)
end

if isfield(opts,'xLimits')
    xlim(opts.xLimits)
else
    xlim([-0.2 nBars+0.2]+0.5)
end

plot(xlim,[opts.baseLine opts.baseLine],'-',...
    'color',[0 0 0],'linewidth',2)

if isfield(opts,'xTicks')
    set(gca,'XTick',opts.xTicks)
else
    set(gca,'XTick',[])
end
if isfield(opts,'XTickLabel')
    set(gca, 'XTickLabel',opts.XTickLabel)
end
if isfield(opts,'yTicks')
    set(gca,'YTick',opts.yTicks)
else
    set(gca,'yTick',[])
end

set(gca,'linewidth',2)
set(gca,'fontweight','bold')
set(gca,'fontsize',15);

if isfield(opts,'yLimits')
    minY = opts.yLimits(1);
    maxY = opts.yLimits(2);
    minY = minY + 0.05*abs(minY);
    ylim([minY maxY]) 
end

if isfield(opts,'save')
    if (opts.save==1) && isfield(opts,'fileName') && isfield(opts,'savePath')
        cPath = pwd;
        cd(opts.savePath)
        addpath(cPath)
        addpath([cPath '/Plotting/'])
        
        try
            fN = strcat(opts.fileName,'.svg');
            %plot2svg(fN,gcf)
            print(gcf,'-dsvg',fN)
            
            eval(['!' inkscapePath ' -z ' fN ' --export-pdf=' opts.fileName '.pdf'])
            eval(['!rm ' fN])
        catch
            keyboard
        end
        cd(cPath)
    end
end

