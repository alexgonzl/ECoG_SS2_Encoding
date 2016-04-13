
function PCAtrialDecomPlots(opts)

% load data and PCA results
fileName = ['allMBAnalysis' opts.lock 'sublogPowernonLPCch'];
load([opts.dataPath opts.lock '/' fileName '.mat'])
fileName = ['PCATrialDecomp-MBAnalysis' opts.lock 'sublogPowernonLPCch'];
load([opts.dataPath opts.lock '/' fileName '.mat'])
pcadat = out; clear out;


    %% plot 1. component correlation histograms
    
    p1Name = ['comp-RTcorr_' opts.lock '_rThr' strrep(num2str(opts.rThr),'.','p')];
    x=[];
    x{1} = pcadat.CorrStudyRTs(:);
    x{2} = pcadat.CorrTestRTs(:);
    titles = {'StudyRTs','TestRTs'};
    figure(1); clf;
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'position',[100,100,600,300]);
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
                if corrBin(bb)>=opts.rThr
                    b(1).FaceColor=[0.9 0.78 0.2];
                elseif corrBin(bb)<=-opts.rThr
                    b(1).FaceColor=[0.2 0.6 0.8];
                end
            end
        end
        set(gca,'fontsize',20,'xTick',[-0.6 -opts.rThr 0 opts.rThr 0.6])
        xlabel('Correlation (r)')
        ylabel('nComps')
        title(titles{tt})
    end
    print(gcf,'-dpdf', [opts.savePath p1Name ])
    
    %% plot2 pcs correlations between study and test
    p2Name = ['study-test_corrs' opts.lock '_rThr' strrep(num2str(opts.rThr),'.','p')];
    
    figure(2); clf;
    x = pcadat.CorrStudyRTs(:);
    y = pcadat.CorrTestRTs(:);
    
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'position',[100,100,500,400]); hold on;
    
    s=scatter(x,y);
    s.MarkerFaceAlpha=0.6;
    s.MarkerEdgeAlpha=0.2;
    s.MarkerEdgeColor = [0.45 0.75 0.9];
    s.MarkerFaceColor = [0.45 0.75 0.9];
    set(gca,'fontsize',20,'xTick',[-0.6 -0.2 0 0.2 0.6],'yTick',[-0.6 -0.2 0 0.2 0.6])
    pfit = polyfit(x,y,1);
    p=plot([-0.5 0.5],[pfit(1)*-0.5 pfit(2)*0.5]+pfit(2));
    p.Color = [0.3 0.3 0.35];
    p.LineWidth=3;
    
    xlabel('Study (r)')
    ylabel('Test (r)')
    xlim([-0.7 0.7])
    print(gcf,'-dpdf', [opts.savePath p2Name ])

%% plot 3 and 4 spectrograms of top components
% component correlations.
x=[];
x{1} = pcadat.CorrStudyRTs;
x{2} = pcadat.CorrTestRTs;

chanIDs=cell(2);
compIDs=cell(2);

%-> 1st row positive correlations for study (col1) and test (col2)
%-> 2nd row negative correlations for study and test
for tt = 1:2
    [chanIDs{1,tt},compIDs{1,tt}]=find(x{tt}>=opts.rThr);
    [chanIDs{2,tt},compIDs{2,tt}]=find(x{tt}<=-opts.rThr);
end

nBins   = sum(data.AnalysisBins);
nBands  = data.nBands;
Xt      = zeros(2,2,data.nBands,nBins);
Xp      = zeros(2,2,data.nBands,nBins);

    for ii = 1:2
        for jj = 1:2
            nComps = numel(chanIDs{ii,jj});
            X    = zeros(nComps,data.nBands*nBins);
            for cc = 1:nComps
                X(cc,:) = pcadat.Projections(chanIDs{ii,jj}(cc),:,compIDs{ii,jj}(cc));
            end
            [~,p,~,t] = ttest(X);
            Xt(ii,jj,:,:) = reshape(t.tstat,[nBands,nBins]);
            Xt2 = squeeze(Xt(ii,jj,:,:));
            
            Xp(ii,jj,:,:) = reshape(p,[nBands,nBins]);
            % threshold activity map after FDR correction
            Yt=reshape(mafdr(p),[nBands nBins]);
            Xt2(Yt>opts.pThr)=0;
            Xt(ii,jj,:,:) = Xt2;
            
        end
    end

    titles =[];
titles{1,1} = ['+StudyRTs-PCs n=' num2str(numel(chanIDs{1,1}))];
titles{1,2} = ['+TestRTs-PCs n=' num2str(numel(chanIDs{1,2}))];
titles{2,1} = ['-StudyRTs-PCs n=' num2str(numel(chanIDs{2,1}))];
titles{2,2} = ['-TestRTs-PCs n=' num2str(numel(chanIDs{2,2}))];

freqs = 1:nBands;
time  = mean(data.Bins(data.AnalysisBins,:),2);

figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1200 1200],'paperposition',[0 0 1200 1200])
set(gcf,'position',[100,100,1000,1000]); % 100pt all around margin

% set axes for each subplot
a = cell(2);
a{1,1} = axes('units','points','position',[50 500 400 300]);      % top left
a{1,2} = axes('units','points','position',[550 500 400 300]);     % top right
a{2,1} = axes('units','points','position',[50 100 400 300]);      % bottom left
a{2,2} = axes('units','points','position',[550 100 400 300]);     % bottom right

% colormap settings
clev = [-4:-1 1:4];
y = brewermap(numel(clev)-1,'*RdBu');
y2 = brewermap(1000,'*RdBu');

% axis tick labels
switch opts.lock
    case 'stim'
        xtick      = [0 0.4 0.8 1.2];
        xticklabel = {'stim','0.4','0.8','1.2'};
    case 'RT'
        xtick      = [-0.8 -0.4 0 0.4];
        xticklabel = {'-0.8','-0.4','resp','0.4'};
    case  'preStim2'
        xtick      = [-.8 -0.4 0 0.4 0.8 1.2];
        xticklabel = {'-0.8','-0.4','stim','0.4','0.8','1.2'};
end
yticksLabels   = {'\delta','\theta','\alpha','\beta','l-\gamma','h-\gamma'};

for ii = 1:2
    for jj=1:2
        axes(a{ii,jj});
        h=contourfcmap(time,freqs,squeeze(Xt(ii,jj,:,:)),clev,y,'lo',y2(1,:),'hi', y2(end,:),'method','calccontour');
        for i = 1:numel(h.l)
            h.l(i).Color='none';
        end
        hold on;
        plot([0 0],[1 6],'k','linewidth',4)
        set(gca,'fontsize',20,'ytick',freqs,'ytickLabel',yticksLabels)
        set(gca,'XTick',xtick,'xtickLabel',xticklabel)
        title(titles{ii,jj})
    end
end
p3Name = ['sigCorrComps' opts.lock '_rThr' strrep(num2str(opts.rThr),'.','p')];
print(gcf,'-dpdf', [opts.savePath p3Name])
print(gcf,'-dtiff','-r300', [opts.savePath p3Name])
