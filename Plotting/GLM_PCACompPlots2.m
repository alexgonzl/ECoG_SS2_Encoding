
function GLM_PCACompPlots2(opts)

% load data and PCA results
fileName = ['allMBAnalysis' opts.lock 'sublogPowernonLPCch'];
load([opts.dataPath opts.lock '/' fileName '.mat'])
fileName = ['PCATrialDecomp-MBAnalysis2_Kmeans' opts.lock 'sublogPowernonLPCch'];
load([opts.dataPath opts.lock '/' fileName '.mat'])
pcadat = out; clear out;
% load electrode locations
load('~/Google Drive/Research/ECoG_SS2e/data_results/Renderings/electrodeLocs.mat')

rThr    = opts.rThr;
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

%% plot 1. Tval histograms
if opts.plot1
    p1Name = ['SelCompsCorr' opts.lock '_rThr' strrep(num2str(rThr),'.','p')];
    x=[];
    x{1} = pcadat.CorrStudyRTs(:);
    x{2} = pcadat.CorrTestRTs(:);
    
    titles = {'StudyRTs','TestRTs'};
    figure(1); clf;
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
    set(gcf,'position',[100,100,600,400]);
    a = axes('position',[0.2 0.2 0.75 0.4]);
    
    for tt=1:2
        cla; axes(a);
        ylim([-1 180]);
        xlim([-0.7 0.7])
        [barHeight,corrBin]=hist(x{tt},-0.7:0.025:0.7);
        nCorrBins = numel(corrBin);
        
        hold on;
        for bb = 1:nCorrBins
            if barHeight(bb)>0
                b=bar(corrBin(bb),barHeight(bb),0.02,...
                    'edgecolor','none','faceColor',[0.85 0.85 0.9]);
                if corrBin(bb)>=(rThr-0.01)
                    b(1).FaceColor=[0.9 0.78 0.2];
                elseif corrBin(bb)<=(-rThr+0.01)
                    b(1).FaceColor=[0.2 0.6 0.8];
                end
                b.ShowBaseLine='off';
            end
        end
        ylabel(' nComps ')
        set(gca,'fontsize',20,'xTick',[-0.5 -rThr 0 rThr 0.5],'linewidth',2)
        xlabel(' Corr (r) ')
        %title(titles{tt})
        print(gcf,'-dpdf', [opts.savePath titles{tt} p1Name ])
    end
    
end

%% plot2 pcs correlations between study and test
if opts.plot2
    p2Name = ['SelCompsCorr_studytest' opts.lock '_rThr' strrep(num2str(rThr),'.','p')];
    x = pcadat.CorrStudyRTs(:);
    y = pcadat.CorrTestRTs(:);
    roisComp = repmat(data.ROIid,[1,pcadat.nComps]);
    roisComp = roisComp(:);
    
    figure(2); clf;
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'paperUnits','points','papersize',[400 400],'paperposition',[0 0 400 400])
    set(gcf,'position',[100,100,400,400]);
    hold on;
    for rr = 1:3
        xx = x(roisComp==rr);
        yy = y(roisComp==rr);
        s=scatter(xx,yy);
        s.MarkerFaceAlpha=0.6;
        s.MarkerFaceColor = ROIcolors(rr,:);
        s.MarkerEdgeColor='none';
    end
    
    xlim([-0.7 0.7])
    ylim([-0.7 0.7])
    set(gca,'fontsize',20,'LineWidth',2)
    xlabel(' StudyRTs (r) ')
    ylabel(' TestRTs (r) ')
    axis square
    print(gcf,'-dpdf', [opts.savePath p2Name ])
end
%% plot 3: Spectrogram of study and test kmeans
if opts.plot3
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
    
    % Get data to plot
    K(1) = pcadat.StudyKmeanOptK;
    K(2) = pcadat.TestKmeanOptK;
    strs = {'StudyCompKmeans','TestCompKmeans'};
    T = cell(2,max(K));
    ROI_CL_Dist = nan(2,max(K),nROIs);
    ids{1} = pcadat.StudySelComps.CompIDs;
    ids{2} = pcadat.TestSelComps.CompIDs;
    for ii =1:2
        for kk = 1:K(ii)
            x = reshape(pcadat.(strs{ii}).T(kk,:), [nBands nBins]);
            p =  pcadat.(strs{ii}).P(kk,:);
            p = reshape(mafdr(p),[nBands,nBins]);
            x(p>opts.pThr) = 0;
            T{ii,kk} = x;
            
            ROI_CL_Dist(ii,kk,:)=histc(rois(ids{ii}(pcadat.(strs{ii}).IDX==kk)),1:nROIs);
        end
    end
    
    figure(4); clf;
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
    set(gcf,'position',[100,100,600,400]); % 100pt all around margin
    a(1) = axes('position',[0.2 0.2 0.6 0.6]);
    a(2) = axes('position',[0.82 0.2 0.1 0.3]);
    for ii = 1:2
        ybarlim=max(nansum(ROI_CL_Dist(ii,:,:),3));
        for kk = 1:K(ii)
            axes(a(1)); cla;
            h=contourfcmap(time,freqs, T{ii,kk},clev,y,'lo',rcolmap(1,:),'hi', rcolmap(end,:),'method','calccontour');
            for i = 1:numel(h.l)
                h.l(i).Color='none';
            end
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            set(gca,'fontsize',20,'ytick',freqs,'ytickLabel',yticksLabels)
            set(gca,'XTick',xtick,'xtickLabel',xticklabel)
            
            nr = squeeze(ROI_CL_Dist(ii,kk,:));
            axes(a(2)); cla
            singleStackBar(nr,ROIcolors,ybarlim,gca);
            
            p3Name = [strs{ii} '_CL' num2str(kk) opts.lock '_pThr' strrep(num2str(opts.pThr),'.','p')];
            print(gcf,'-dtiff','-r300', [opts.savePath p3Name])
        end
    end
end

%% plot 4: Spectrogram of study and test kmeans by region.
if opts.plot4
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
    
    strs = {'StudyKmeanCompROIs','TestKmeanCompROIs'};
    figure(5); clf;
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
    set(gcf,'position',[100,100,600,400]); % 100pt all around margin
    a   = axes('position',[0.15 0.15 0.70 0.5]);
    
    K = nan(2,nROIs);
    for ii = 1:2
        K(ii,:)=pcadat.(strs{ii}).Ko;
        for rr = 1:nROIs
            Kr = K(ii,rr);
            for kk =1:Kr
                x = reshape(pcadat.(strs{ii}).T{rr}(kk,:),[nBands,nBins]);
                p = pcadat.(strs{ii}).P{rr}(kk,:);
                p = reshape(mafdr(p),[nBands,nBins]);
                x(p>opts.pThr)  = 0;
                
                axes(a(1)); cla;
                h=contourfcmap(time,freqs, x,clev,y,'lo',rcolmap(1,:),'hi', rcolmap(end,:),'method','calccontour');
                for i = 1:numel(h.l)
                    h.l(i).Color='none';
                end
                hold on;
                plot([0 0],[1 6],'k','linewidth',4)
                set(gca,'fontsize',20,'ytick',freqs,'ytickLabel',yticksLabels)
                set(gca,'XTick',xtick,'xtickLabel',xticklabel)
                pName = [strs{ii} '_' data.ROIs{rr} '_CL' num2str(kk) opts.lock '_pThr' strrep(num2str(opts.pThr),'.','p')];
                if ~exist([opts.savePath '/ROI_CL/'],'dir')
                    mkdir([opts.savePath '/ROI_CL/'])
                end
                print(gcf,'-dtiff','-r300', [opts.savePath '/ROI_CL/' pName])
            end
        end
    end
end

%% plot 5: Spectrograms of PC of selected
if opts.plot5
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
    nCompsToPlot = 3;
    strs = {'StudyPCASelComps','TestPCASelComps'};
    figure(5); clf;
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
    set(gcf,'position',[100,100,600,400]); % 100pt all around margin
    a   = axes('position',[0.15 0.15 0.70 0.5]);
    
    for ii=1:2
        for jj = 1:nCompsToPlot
            x = reshape(pcadat.(strs{ii}).S(:,jj),[nBands nBins]);
            v = round(pcadat.(strs{ii}).E(jj));
            axes(a); cla;
            clev2 = quantile(x(:),(clev+5)/10);
            h=contourfcmap(time,freqs, x,clev2,y,'lo',rcolmap(1,:),'hi', rcolmap(end,:),'method','calccontour');
            for i = 1:numel(h.l)
                h.l(i).Color='none';
            end
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            set(gca,'fontsize',20,'ytick',freqs,'ytickLabel',yticksLabels)
            set(gca,'XTick',xtick,'xtickLabel',xticklabel)
            title(sprintf('VarExp %d',v))
            pName = [strs{ii} '_Comp' num2str(jj) opts.lock];            
            print(gcf,'-dtiff','-r300', [opts.savePath pName])
            
        end
    end
end
%% spectrograms of PC of componenets by region
if opts.plot6
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
    nCompsToPlot = 3;
    strs = {'StudyPCASelCompsROIs','TestPCASelCompsROIs'};
    figure(5); clf;
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
    set(gcf,'position',[100,100,600,400]); % 100pt all around margin
    a   = axes('position',[0.15 0.15 0.70 0.5]);
    for ii=1:2
        for rr = 1:3
            for jj = 1:nCompsToPlot
                x = reshape(pcadat.(strs{ii}).S{rr}(:,jj),[nBands nBins]);
                v = round(pcadat.(strs{ii}).E{rr}(jj));
                axes(a); cla;
                clev2 = quantile(x(:),(clev+5)/10);
                h=contourfcmap(time,freqs, x,clev2,y,'lo',rcolmap(1,:),'hi', rcolmap(end,:),'method','calccontour');
                for i = 1:numel(h.l)
                    h.l(i).Color='none';
                end
                hold on;
                plot([0 0],[1 6],'k','linewidth',4)
                set(gca,'fontsize',20,'ytick',freqs,'ytickLabel',yticksLabels)
                set(gca,'XTick',xtick,'xtickLabel',xticklabel)
                title(sprintf('VarExp %d',v))
                pName = [strs{ii} '_' data.ROIs{rr} '_Comp' num2str(jj) opts.lock];
                if ~exist([opts.savePath '/ROI_Comps/'],'dir')
                    mkdir([opts.savePath '/ROI_Comps/'])
                end
                print(gcf,'-dtiff','-r300', [opts.savePath '/ROI_Comps/' pName])
            end
        end
    end
    
end
end
%% aux funcs
function han=singleStackBar(x,cols,ybarlim,han)

x = x(:);
nCat = numel(x);
if nargin<3
    han = axes('position',[0.1 0.1 0.2 0.8]);
end

axes(han); cla; hold on;
xx = cumsum([x]);
for ii = nCat:-1:1
    b=bar(1,xx(ii));
    b.FaceColor =cols(ii,:);
    b.EdgeColor='none';
    b.ShowBaseLine='off';
end
ylim([0 ybarlim])
set(gca,'xtick',[])
axis off
end

