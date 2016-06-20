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
ha=figure(); clf;hold on;
set(gcf,'paperpositionmode','auto','position',[100 100 opts.aspectRatio])

for ba = 1:nBars
    
    X = M{ba};
    n = size(X,1);
    [m(ba) se(ba)] = grpstats(X,ones(n,1),{'mean','sem'});
    
    bar(opts.xPositions(ba),m(ba),0.8,'FaceColor',opts.colors(ba,:),'edgeColor', 'none', 'basevalue',opts.baseLine,'ShowBaseLine','off')
    plot([1 1]*opts.xPositions(ba), [-se(ba) se(ba)]+m(ba), 'color',[0 0 0],'linewidth',4)
    
end
ma = max(m+se);

plot(xlim,[opts.baseLine opts.baseLine],'-',...
    'color',[0 0 0],'linewidth',2)

if isfield(opts,'xTick')
    set(gca,'xTick',opts.xTick)
else
    set(gca,'xTick',[])
end

if isfield(opts,'xTickLabel')
    set(gca,'xTickLabel',opts.xTickLabel)
else
    set(gca,'xTickLabel',[])
end

if isfield(opts,'yTicks')
    set(gca,'YTick',opts.yTicks)
else
    set(gca,'yTick',[])
end

if isfield(opts,'yLimits')
    ylim(opts.yLimits)
else        
    ylim([0 ma*1.25])    
end

if isfield(opts,'xLimits')
    xlim(opts.xLimits)
else
    xlim([0 nBars]+0.5)
end

if isfield(opts,'yLabel')
    ylabel(opts.yLabel)
end
if isfield(opts, 'sigMarks')
    nMarks = size(opts.sigMarks,1);
    if isfield(opts,'sigLevel')
        marks=cell(nMarks,1);
        for ii=1:nMarks,
            if opts.sigLevel(ii)<0.001
                marks{ii}='***';
            elseif opts.sigLevel(ii)<0.01
                marks{ii}='**';
            elseif opts.sigLevel(ii)<0.05
                marks{ii}='*';
            elseif opts.sigLevel(ii)<0.1
                marks{ii}='~';
            end
        end
    else
        marks=cell(nMarks,1);
        for ii=1:nMarks, marks{ii}='*'; end
    end
    for ii = 1:nMarks
        ypos = [1 1]*ma*(1+.18*ii);
        xpos = opts.sigMarks(ii,:);
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
