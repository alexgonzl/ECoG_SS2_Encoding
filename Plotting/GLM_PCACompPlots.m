
function GLM_PCACompPlots(opts)

% load data and PCA results
fileName = ['allMBAnalysis' opts.lock 'sublogPowernonLPCch'];
load([opts.dataPath opts.lock '/' fileName '.mat'])
fileName = ['PCATrialDecomp-MBAnalysis' opts.lock 'sublogPowernonLPCch'];
load([opts.dataPath opts.lock '/' fileName '.mat'])
pcadat = out; clear out;
% load electrode locations
load('~/Google Drive/Research/ECoG_SS2e/data_results/Renderings/electrodeLocs.mat')

tThr    = pcadat.GLMsCompsThr;
freqs   = 1:data.nBands;
time    = mean(data.Bins(data.AnalysisBins,:),2);
nFeat   = pcadat.nFeat;
nChans  = pcadat.nChans;
rois    = data.ROIid;
nROIs   = numel(unique(rois));
% colormap settings
clev = [-4:-1 1:4];
y = brewermap(numel(clev)-1,'*RdBu');
rcolmap = brewermap(1000,'*RdBu');
nBins   = sum(data.AnalysisBins);
ROIcolors = [240 35 17; 2 93 140;122 179 23]/255;  

%% plot 1. Tval histograms
if opts.plot1
    p1Name = ['GLM_PCA-CompTs_' opts.lock '_tThr' strrep(num2str(tThr),'.','p')];
    x=[];
    x{1} = pcadat.StudyGLMSChanCompTVal(:);
    x{2} = pcadat.TestGLMSChanCompTVal(:);
    
    titles = {'StudyRTs','TestRTs'};
    figure(1); clf;
    set(gcf,'paperpositionmode','auto','color','white')    
    set(gcf,'paperUnits','points','papersize',[800 500],'paperposition',[0 0 800 600])
    set(gcf,'position',[100,100,600,400]);
    a{1} = axes('position',[0.1 0.2 0.4 0.4]);
    a{2} = axes('position',[0.55 0.2 0.4 0.4]);
    
    for tt=1:2
        axes(a{tt})
        ylim([0 80]);
        xlim([-4 4])
        [barHeight,corrBin]=hist(x{tt},-4:0.1:4);
        nCorrBins = numel(corrBin);
        
        hold on;
        for bb = 1:nCorrBins
            if barHeight(bb)>0
                b=bar(corrBin(bb),barHeight(bb),5/nCorrBins,...
                    'edgecolor','none','faceColor',[0.85 0.85 0.9]);
                if corrBin(bb)>=tThr
                    b(1).FaceColor=[0.9 0.78 0.2];
                elseif corrBin(bb)<=-tThr
                    b(1).FaceColor=[0.2 0.6 0.8];
                end
            end
        end
        if tt==1
            ylabel(' nComps ')
        else
            set(gca,'YTickLabel',[])
        end
        set(gca,'fontsize',20,'xTick',[-3 -tThr 0 tThr 3],'linewidth',2)
        xlabel(' T ')
        title(titles{tt})
    end
    print(gcf,'-dpdf', [opts.savePath p1Name ])
end

%% plot2 pcs correlations between study and test
if opts.plot2
    p2Name = ['GLM_PCA-Comp_studytestCorr' opts.lock '_tThr' strrep(num2str(tThr),'.','p')];
    
    figure(2); clf;
    x = pcadat.StudyGLMSChanCompTVal(:);
    y = pcadat.TestGLMSChanCompTVal(:);
    
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'position',[100,100,500,400]); hold on;
    %
    s=scatter(x,y);
    xlim([-4.5 4.5])
    ylim([-4.5 4.5])
    s.MarkerFaceAlpha=0.6;
    s.MarkerEdgeAlpha=0.2;
    s.MarkerEdgeColor = [0.45 0.75 0.9];
    s.MarkerFaceColor = [0.45 0.75 0.9];
    %set(gca,'fontsize',20,'xTick',[0 opts.rThr 0.6],'yTick',[ 0 opts.rThr 0.6])
    set(gca,'fontsize',20)
    pfit = polyfit(x,y,1);
    p=plot([-3 3],[pfit(1)*-3 pfit(1)*3]+pfit(2));
    p.Color = [0.3 0.3 0.35];
    p.LineWidth=3;
    %
    xlabel(' Study (T) ')
    ylabel(' Test (T) ')
    
    print(gcf,'-dpdf', [opts.savePath p2Name ])
end
%% plot 3: Spectrogram of positive and negative components
if opts.plot3
    ChCompIDs = cell(2,2); % rows: pos/neg col: study/test
    ChCompIDs{1,1}= pcadat.StudyGLMsCompKmeans.PosCompIDs;
    ChCompIDs{1,2}= pcadat.TestGLMsCompKmeans.PosCompIDs;
    ChCompIDs{2,1}= pcadat.StudyGLMsCompKmeans.NegCompIDs;
    ChCompIDs{2,2}= pcadat.TestGLMsCompKmeans.NegCompIDs;
    
    tThr    = pcadat.GLMsCompsThr;
    nBands  = data.nBands;
    
    X = cell(2,2);% rows: pos/neg col: study/test
    for ii = 1:2
        for jj=1:2
            nComps = size(ChCompIDs{ii,jj},1);
            Y = zeros(nComps,nFeat);
            for kk = 1:nComps
                ch = ChCompIDs{ii,jj}(kk,1);
                co = ChCompIDs{ii,jj}(kk,2);
                Y(kk,:) = pcadat.Projections(ch,:,co);
            end
            [~,p,~,t] = ttest(Y);
            p = reshape(mafdr(p),[nBands,nBins]);
            x = reshape(t.tstat,[nBands,nBins]);
            x(p>opts.pThr) = 0;
            X{ii,jj} = x;
        end
    end
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
    
    
    figure(4); clf;
    set(gcf,'paperpositionmode','auto','color','white')
    set(gcf,'paperUnits','points','papersize',[1200 1000],'paperposition',[0 0 1200 1000])
    set(gcf,'position',[100,100,1000,800]); % 100pt all around margin
    
    % set axes for each subplot
    a = cell(2,2);
    a{1,1} = axes('units','points','position',[50 400 400 200]);      % mid left
    a{1,2} = axes('units','points','position',[550 400 400 200]);     % mid right
    a{2,1} = axes('units','points','position',[50 100 400 200]);      % bottom left
    a{2,2} = axes('units','points','position',[550 100 400 200]);     % bottom right
    
    for ii = 1:2 % clusters
        for jj=1:2 % pos/neg
            axes(a{ii,jj});
            h=contourfcmap(time,freqs, X{ii,jj},clev,y,'lo',rcolmap(1,:),'hi', rcolmap(end,:),'method','calccontour');
            for i = 1:numel(h.l)
                h.l(i).Color='none';
            end
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            set(gca,'fontsize',20,'ytick',freqs,'ytickLabel',yticksLabels)
            set(gca,'XTick',xtick,'xtickLabel',xticklabel)
            %title(titles{kk,jj,ii})
        end
    end
    p3Name = ['CompsTVals' opts.lock strrep(num2str(opts.pThr),'.','p')];
    print(gcf,'-dtiff','-r300', [opts.savePath p3Name])
end

%% plot 4 spectrograms of top components
% component correlations.
if opts.plot4
    K = pcadat.GLMsCompsKMeans;
    tThr = pcadat.GLMsCompsThr;
    nBins   = sum(data.AnalysisBins);
    nBands  = data.nBands;
    
    X = cell(K,2,2);
    str1 = {'StudyGLMsCompKmeans','TestGLMsCompKmeans'};
    str2 = {'PT','PP';'NT','NP'};
    for kk = 1:K
        for ii = 1:2
            for jj=1:2
            x = reshape(pcadat.(str1{ii}).(str2{jj,1})(kk,:),[nBands,nBins]);
            p = reshape(mafdr(pcadat.(str1{ii}).(str2{jj,2})(kk,:)),[nBands,nBins]);
            x(p>opts.pThr) = 0;
            X{kk,ii,jj} = x;
            end
        end
    end
    
    titles  = cell(K,2,2);
    Ns      = zeros(K,2,2);
    str1 = {'StudyGLMsCompKmeans','TestGLMsCompKmeans'};
    str2 = {'PIDX','NIDX'};
    str3 = {'+StudyPCs','-StudyPCs'; '+TestPCs','-TestPCs'};
    for ii =1:2 % study/test
        for jj=1:2 % pos/neg
            Ns(:,jj,ii)     = histc(pcadat.(str1{ii}).(str2{jj}),1:K);
            for kk=1:K % clusters
                titles{kk,jj,ii} = ['K=',num2str(kk) , ', ' str3{ii,jj}, ' n=' num2str(Ns(kk,jj,ii))];
            end
        end
    end
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
    
    for ii = 1:2 % study /test
        figure(ii+2); clf;
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'paperUnits','points','papersize',[1200 1200],'paperposition',[0 0 1200 1200])
        set(gcf,'position',[100,100,1000,1000]); % 100pt all around margin
        
        % set axes for each subplot
        a = cell(3,2);
        a{1,1} = axes('units','points','position',[50 700 400 200]);      % mid left
        a{1,2} = axes('units','points','position',[550 700 400 200]);     % mid right
        a{2,1} = axes('units','points','position',[50 400 400 200]);      % mid left
        a{2,2} = axes('units','points','position',[550 400 400 200]);     % mid right
        a{3,1} = axes('units','points','position',[50 100 400 200]);      % bottom left
        a{3,2} = axes('units','points','position',[550 100 400 200]);     % bottom right
        
        for kk = 1:K % clusters
            for jj=1:2 % pos/neg
                axes(a{kk,jj});
                h=contourfcmap(time,freqs, X{kk,jj,ii},clev,y,'lo',rcolmap(1,:),'hi', rcolmap(end,:),'method','calccontour');
                for i = 1:numel(h.l)
                    h.l(i).Color='none';
                end
                hold on;
                plot([0 0],[1 6],'k','linewidth',4)
                set(gca,'fontsize',20,'ytick',freqs,'ytickLabel',yticksLabels)
                set(gca,'XTick',xtick,'xtickLabel',xticklabel)
                title(titles{kk,jj,ii})
            end
        end
        p4Name = [str1{ii} 'Clusters' opts.lock strrep(num2str(opts.pThr),'.','p')];
        print(gcf,'-dtiff','-r300', [opts.savePath p4Name])
    end
end

%% plot 5: renderings of top components

if opts.plot5

    % level indicates the # of componnents above thr by channel
    % rows: pos/neg col:study/test
    strs = {'PosStudyChRend','PosTestChRend';'NegStudyChRend','NegTestChRend'};
    Weights = [];
    Weights{1,1}=sum(pcadat.StudyGLMSChanCompTVal>tThr,2);
    Weights{1,2}=sum(pcadat.TestGLMSChanCompTVal>tThr,2);
    Weights{2,1}=sum(pcadat.StudyGLMSChanCompTVal<-tThr,2);
    Weights{2,2}=sum(pcadat.TestGLMSChanCompTVal<-tThr,2);
    
    levels   = 1:8;
    nLevels = numel(levels);
    for ii = 1:2
        if ii==1
            CM =  brewermap(nLevels,'YlOrRd');
        else
            CM =  brewermap(nLevels,'YlGnBu');
        end
        for jj = 1:2
            w = Weights{ii,jj};
            w(w==0) = nan;
            han=renderChanWeights(elecLocs,w+1,CM);
            print(han,'-dtiff','-r300', [opts.savePath strs{ii,jj}])
        end
    end
    
end
%% plot 6: same as plot 4 set (spectrograms), with relative contributions by region
if opts.plot6
    K = pcadat.GLMsCompsKMeans;
    nBins   = sum(data.AnalysisBins);
    nBands  = data.nBands;
    
    X = cell(K,2,2);
    str1 = {'StudyGLMsCompKmeans','TestGLMsCompKmeans'};
    str2 = {'PT','PP';'NT','NP'};
    for kk = 1:K
        for ii = 1:2
            for jj=1:2
            x = reshape(pcadat.(str1{ii}).(str2{jj,1})(kk,:),[nBands,nBins]);
            p = reshape(mafdr(pcadat.(str1{ii}).(str2{jj,2})(kk,:)),[nBands,nBins]);
            x(p>opts.pThr) = 0;
            X{kk,ii,jj} = x;
            end
        end
    end
    
    titles  = cxell(K,2,2);
    Ns      = zeros(K,2,2);
    str1 = {'StudyGLMsCompKmeans','TestGLMsCompKmeans'};
    str2 = {'PIDX','NIDX'};
    str3 = {'+StudyPCs','-StudyPCs'; '+TestPCs','-TestPCs'};
    for ii =1:2 % study/test
        for jj=1:2 % pos/neg
            Ns(:,jj,ii)     = histc(pcadat.(str1{ii}).(str2{jj}),1:K);
            for kk=1:K % clusters
                titles{kk,jj,ii} = ['K=',num2str(kk) , ', ' str3{ii,jj}, ' n=' num2str(Ns(kk,jj,ii))];
            end
        end
    end
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
    
    for ii = 1:2 % study /test
        
        figure(ii+2); clf;
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'paperUnits','points','papersize',[1300 1200],'paperposition',[0 0 1200 1200])
        set(gcf,'position',[100,100,1100,1000]); % 100pt all around margin
        
        % set axes for each subplot
        a = cell(3,2,2);
        a{1,1,1} = axes('units','points','position',[50 700 400 200]);      % top left
        a{1,1,2} = axes('units','points','position',[460 700 50 100]);      % top leftb
        a{1,2,1} = axes('units','points','position',[600 700 400 200]);     % top right
        a{1,2,2} = axes('units','points','position',[1010 700 50 100]);     % top rightb
        a{2,1,1} = axes('units','points','position',[50 400 400 200]);      % mid left
        a{2,1,2} = axes('units','points','position',[460 400 50 100]);      % mid leftb        
        a{2,2,1} = axes('units','points','position',[600 400 400 200]);     % mid right
        a{2,2,2} = axes('units','points','position',[1010 400 50 100]);     % mid rightb
        a{3,1,1} = axes('units','points','position',[50 100 400 200]);      % bottom left
        a{3,1,2} = axes('units','points','position',[460 100 50 100]);      % bottom leftb        
        a{3,2,1} = axes('units','points','position',[600 100 400 200]);     % bottom right
        a{3,2,2} = axes('units','points','position',[1010 100 50 100]);     % bottom rightb
        
        ybarlim=max(max(Ns(:,:,ii)));
        for kk = 1:K % clusters
            for jj=1:2 % pos/neg
                
                axes(a{kk,jj,1});
                h=contourfcmap(time,freqs, X{kk,jj,ii},clev,y,'lo',rcolmap(1,:),'hi', rcolmap(end,:),'method','calccontour');
                for i = 1:numel(h.l)
                    h.l(i).Color='none';
                end
                hold on;
                plot([0 0],[1 6],'k','linewidth',4)
                set(gca,'fontsize',20,'ytick',freqs,'ytickLabel',yticksLabels)
                set(gca,'XTick',xtick,'xtickLabel',xticklabel)
                title(titles{kk,jj,ii})
                                
                nr = pcadat.GLMsThrComps.nCompsByCLROI(:,kk,ii,jj) ;
                singleStackBar(nr,ROIcolors,ybarlim,a{kk,jj,2});
            end
        end
        p6Name = [str1{ii} 'Clusters2' opts.lock strrep(num2str(opts.pThr),'.','p') '_tThr' strrep(num2str(tThr),'.','p')];
        print(gcf,'-dtiff','-r300', [opts.savePath p6Name])
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

