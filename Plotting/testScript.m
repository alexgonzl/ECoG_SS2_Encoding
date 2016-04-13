
function PCAtrialDecomPlots(data,opts)

%%
% plot 1. component correlation histograms
x=[];
x{1} = data.CorrStudyRTs(:);
x{2} = data.CorrTestRTs(:);
titles = {'StudyRTs','TestRTs'};
figure(1); clf;
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'position',[100,100,1000,400]);
    a{1} = axes('position',[0.1 0.15 0.4 0.7]);
    a{2} = axes('position',[0.55 0.15 0.4 0.7]);
    
for tt=1:2
    axes(a{tt})
    ylim([0 250]);
    xlim([-0.7 0.7])
    [barHeight,corrBin]=hist(x{tt},-1:0.05:1);    
    nCorrBins = numel(corrBin);

    hold on;
    for bb = 1:nCorrBins
        if barHeight(bb)>0
            b=bar(corrBin(bb),barHeight(bb),1/nCorrBins*2,...
                'edgecolor','none','faceColor',[0.85 0.85 0.9]);
            if corrBin(bb)>=0.2
                b(1).FaceColor=[0.9 0.78 0.2];
            elseif corrBin(bb)<=-0.2
                b(1).FaceColor=[0.2 0.6 0.8];
            end
        end
    end
    set(gca,'fontsize',20,'xTick',[-0.6 -0.2 0 0.2 0.6])
    xlabel('Correlation (r)')
    ylabel('nComps')
    title(titles{tt})
end
print(gcf,'-dpdf', [opts.savePath ])

% plot2 