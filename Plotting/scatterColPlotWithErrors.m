function ha = scatterColPlot(M,opts)
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

nBars   = size(M,1);
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

if ~isfield(opts,'markerSize')
    opts.markerSize = 25;
end

ha=figure(); clf;hold on;
set(gcf,'paperpositionmode','auto','position',[100 100 opts.aspectRatio])

for ba = 1:nBars
    
    X = M{ba};
    n = size(X,1);
    [m(ba), ci(ba)] = grpstats(X,ones(n,1),{'mean','ci'});

    xpos = opts.xPositions(ba)+randn(n,1)*0.1;
    s = scatter(xpos,X,opts.markerSize,opts.colors(ba,:),'filled');
    plot([-0.25 0.25]+opts.xPositions(ba), [1 1]*m(ba), 'color',opts.colors(ba,:),'linewidth',4)
    plot([1 1]*opts.xPositions(ba), [-se(ba) se(ba)]+m(ba), 'color',[0.1 0.1 0.1],'linewidth',4)
    
end
ma = cellfun(@max,M); ma = max(ma(:));
mi = cellfun(@min,M); mi = max(mi(:));

plot(xlim,[opts.baseLine opts.baseLine],'-',...
    'color',[0 0 0],'linewidth',2)

if isfield(opts,'xTick')
    set(gca,'xTick',opts.xTick)
else
    set(gca,'xTick',[])
end

if isfield(opts,'grid')
    grid on;
end

if isfield(opts,'xTickLabel')
    set(gca,'xTickLabel',opts.xTickLabel)
else
    set(gca,'xTickLabel',[])
end

if isfield(opts,'yTicks')
    if strcmp(opts.yTicks,'def')
        
    else
        set(gca,'YTick',opts.yTicks)
    end
else
    set(gca,'yTick',[])
end

if isfield(opts,'yLimits')
    ylim(opts.yLimits)
else
    ylim([min(0,mi*1.4) max(0,ma*1.4)])    
end

if isfield(opts,'xLimits')
    xlim(opts.xLimits)
else
    xlim([0 nBars]+0.5)
end

if isfield(opts,'yLabel')
    ylabel(opts.yLabel)
end

if isfield(opts, 'sigUniMarks')
    nMarks = size(opts.sigUniMarks,1);
    if isfield(opts,'sigUniLevel')
        marks=cell(nMarks,1);
        for ii=1:nMarks
            if opts.sigUniLevel(ii)<0.001
                marks{ii}='***';
            elseif opts.sigUniLevel(ii)<0.01
                marks{ii}='**';
            elseif opts.sigUniLevel(ii)<0.05
                marks{ii}='*';
            elseif opts.sigUniLevel(ii)<0.1
                marks{ii}='~';
            end
        end
    else
        marks=cell(nMarks,1);
        for ii=1:nMarks, marks{ii}='*'; end
    end
    for ii = 1:nMarks
        ypos = 1.05*ma;
        xpos = opts.sigUniMarks(ii);     
        if strcmp(marks{ii},'~')            
            tt=text(xpos,ypos*1.05, marks{ii});
        else
            tt=text(xpos,ypos*1.02, marks{ii});
        end
        tt.VerticalAlignment='middle'; tt.HorizontalAlignment='center';
        tt.FontSize=40;
        %plot(mean(xpos),ma*(1.04+.12*ii),marks{ii},'color','k')
    end    
end

if isfield(opts, 'sigBiMarks')
    nMarks = size(opts.sigBiMarks,1);
    if isfield(opts,'sigBiLevel')
        marks=cell(nMarks,1);
        for ii=1:nMarks
            if opts.sigBiLevel(ii)<0.001
                marks{ii}='***';
            elseif opts.sigBiLevel(ii)<0.01
                marks{ii}='**';
            elseif opts.sigBiLevel(ii)<0.05
                marks{ii}='*';
            elseif opts.sigBiLevel(ii)<0.1
                marks{ii}='~';
            end
        end
    else
        marks=cell(nMarks,1);
        for ii=1:nMarks, marks{ii}='*'; end
    end
    for ii = 1:nMarks
        ypos = [1 1]*ma*(1+.11*ii);
        xpos = opts.sigBiMarks(ii,:);
        plot(xpos,ypos,'-k','linewidth',2)
        if strcmp(marks{ii},'~')            
            tt=text(mean(xpos),ypos(1)*1.05, marks{ii});
        else
            tt=text(mean(xpos),ypos(1)*1.02, marks{ii});
        end
        tt.VerticalAlignment='middle'; tt.HorizontalAlignment='center';
        tt.FontSize=25;
        %plot(mean(xpos),ma*(1.04+.12*ii),marks{ii},'color','k')
    end    
end

set(gca,'LineWidth',2,'FontSize',20)
