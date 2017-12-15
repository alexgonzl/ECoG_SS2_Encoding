function ha = scatterColPlot_grp(M,opts)
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

[nGrp, nBars]  = size(M);
N = numel(M);

m       = zeros(nGrp,nBars);
se      = zeros(nGrp,nBars);

if ~isfield(opts,'aspectRatio')
    opts.aspectRatio = [200 600];
end

if ~isfield(opts,'colors')
    for ii=1:nGrp
        opts.colors{ii} = zeros(nBars,3);
    end
end

if ~isfield(opts,'baseLine')
    opts.baseLine =0;
end

if ~isfield(opts,'xPositions')
    opts.xPositions = zeros(nGrp,nBars);
    
    for ii =1:nGrp
        opts.xPositions(ii,:) = (1:nBars) + (ii-1)*(nBars +1);
    end
end

if ~isfield(opts,'markerSize')
    opts.markerSize = 40;
end

for gg = 1:nGrp
    for ba = 1:nBars
        
        X = M{gg,ba};
        n = size(X,1);
        [m(gg,ba), se(gg,ba)] = grpstats(X,ones(n,1),{'mean','sem'});
        
        
        xpos = opts.xPositions(gg,ba)+randn(n,1)*0.1;
        s = scatter(xpos,X,opts.markerSize,opts.colors{gg}(ba,:),'filled');
        s.MarkerFaceAlpha= 0.8;
        plot([-0.25 0.25]+opts.xPositions(gg,ba), [1 1]*m(gg,ba), 'color',[0 0 0],'linewidth',4)
        %plot([1 1]*opts.xPositions(ba), [-se(ba) se(ba)]+m(ba), 'color',[0.1 0.1 0.1],'linewidth',4)
    end
end
ma = cellfun(@max,M); ma = max(ma(:));
mi = cellfun(@min,M); mi = min(mi(:));

% plot(xlim,[opts.baseLine opts.baseLine],'-',...
%     'color',[0 0 0],'linewidth',2)


if isfield(opts,'xTick')
    if strcmp(opts.xTick,'def')
        xx=opts.xPositions';
        set(gca,'xTick',xx(:))
    else
        set(gca,'xTick',opts.xTick)
    end
else
    set(gca,'xTick',[])
end

if isfield(opts,'grid')
    grid on;
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
    xlim([0 opts.xPositions(nGrp,nBars)]+0.5)
end

if isfield(opts,'yLabel')
    ylabel(opts.yLabel)
end

if isfield(opts, 'sigUniMarks')
    [n1,n2] = size(opts.sigUniMarks);
    if isfield(opts,'sigUniLevel')
        marks=cell(n1,n2);
        for ii=1:n1
            for jj = 1:n2
                if opts.sigUniLevel(ii,jj)<0.001
                    marks{ii,jj}='***';
                elseif opts.sigUniLevel(ii,jj)<0.01
                    marks{ii,jj}='**';
                elseif opts.sigUniLevel(ii,jj)<0.05
                    marks{ii,jj}='*';
                elseif opts.sigUniLevel(ii,jj)<0.1
                    marks{ii,jj}='~';
                end
            end
        end
    else
        marks=cell(n1,n2);
        for ii=1:n1
            for jj=1:n2
                if opts.sigUniMarks(ii,jj)
                    marks{ii,jj}='*';
                end
            end
        end
    end
    for ii = 1:n1
        for jj=1:n2
            if opts.sigUniMarks(ii,jj)
                ypos = 1.05*ma;
                xpos = opts.xPositions(ii,jj);
                if strcmp(marks{ii,jj},'~')
                    tt=text(xpos,ypos*1.05, marks{ii,jj});
                else
                    tt=text(xpos,ypos*1.02, marks{ii,jj});
                end
                tt.VerticalAlignment='middle'; tt.HorizontalAlignment='center';
                tt.FontSize=40;
            end
        end
        %plot(mean(xpos),ma*(1.04+.12*ii),marks{ii},'color','k')
    end
end

set(gca,'LineWidth',2,'FontSize',20)


if isfield(opts,'xTickLabel')
    set(gca,'xTickLabel',repmat(opts.xTickLabel{1},1,nGrp))
    xx= mean(opts.xPositions,2);
    ax = gca;
    ax.Position(2)=0.15;
    ax.Position(4)=0.8;
    xLims = xlim;
    
    axes('position',[ax.Position(1) 0, ax.Position(3), 0.12])
    ax=gca;
    xlim(xLims); ylim([0 1]);
    for ii =1:nGrp
        tt=text(xx(ii),0.3,opts.xTickLabel{2}(ii));
        tt.VerticalAlignment='middle'; tt.HorizontalAlignment='center';
        tt.FontSize=25;
    end
    ax.Color='none';
    ax.YColor='none';
    ax.XColor='none';
else
    set(gca,'xTickLabel',[])
end
end
