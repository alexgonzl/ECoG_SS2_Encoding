function SS2e_MBAnalysis(lock,AnalysisNum,opts)
inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';

addpath Analysis/
addpath Plotting/
addpath lib/

%dirPath         = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
dirPath         = ['~/Google Drive/Research/ECoG_SS2e/data_results/'] ;
%fileName        = ['allERSPs' band 'Group' lock 'sublogPowernonLPCCh'];
%fileName        = ['allMBAnalysis' lock 'sublogPowernLPClowEvokedVar'];
fileName        = ['allMBAnalysis_2' lock 'sublogPowernonLPCch'];

load([dirPath lock '/' fileName '.mat'])

% channel groupings/ROIs
nChans = numel(data.ROIid);
groups = zeros(nChans,1);
groups(data.ROIid==1) = 1; %IPS
groups(data.ROIid==2) = 2; %SPL
groups(data.ROIid==3) = 3; %AG


%%
close all

info            =[];
info.groups     = groups;
info.rownames   = cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
info.savePath   = opts.savePath;
% info.savePath   = ['~/Google Drive/Research/ECoG_SS2e' ...
%     '/Plots/SS2e_MB_Analyses/a' num2str(AnalysisNum) '/'];
if ~exist(info.savePath,'dir')
    mkdir(info.savePath);
end

info.groupNames = {'IPS','SPL','AG'};

info.Bins       = data.Bins;
%info.yLimits    = [-1 1.6];
%info.yLimits    = [-4 2];
switch lock
    case 'stim'
    info.xtick      = [0 0.4 0.8 1.2];
    info.toi        = [0 1.5];
    info.alpha      = 0.05/15;
    info.xticklabel = {'stim','0.4','0.8','1.2'};
    case 'RT'
    info.yAxisRightLoc = 1;
    info.xtick      = [-0.8 -0.4 0 0.4];
    info.xticklabel = {'-0.8','-0.4','resp','0.4'};
    case  {'preStim2'}
    info.xtick      = [-.8 -0.4 0 0.4 0.8 1.2];
    info.toi        = [-1 1.5];
    info.alpha      = 0.05/26;
    info.xticklabel = {'-0.8','-0.4','stim','0.4','0.8','1.2'};
    case  {'preStim'}
    %info.xtick      = [-1.2 -.8 -0.4 0 0.4 0.8 1.2];    
    info.toi        = [-1.5 1.5];
    info.alpha      = 0.05/30;
    %info.xticklabel = {'-1.2','-0.8','-0.4','stim','0.4','0.8','1.2'};
    %info.xticklabel = {'-0.8','stim','0.8'};
    info.xtick      = [0 0.5 1]; %%0.5 1];
    info.xticklabel = {'stim','0.5','1'};
end

info.cols            = 'oc';
info.yticksLabels   = {'\delta','\theta','\alpha','\beta','l-\gamma','h-\gamma'};
info.thr            = 1;
info.pThr           = opts.pThr; % corrected.

% Analysis 1

y2 = brewermap(1000,'*RdBu');


switch AnalysisNum
    % Plot ROIs level activity
    case 1
        AnalysisBins = data.Bins(:,1)>-0.51;
        clev = [-5:-1 1:5];
        y = brewermap(numel(clev)-1,'*RdBu');
        Zlims = [-5 5];
        figure(1); clf;
        %        colormap(cm);
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,300]);
        info.Analysis   = [lock '_SpectoROIsTVals'];
        
        freqs = 1:data.nBands; nBands = data.nBands;
        nBins = sum(AnalysisBins);
        time  = mean(data.Bins(AnalysisBins,:),2);
        for rr =1:3
            % select ROI spectrogram
            X = squeeze(data.tBinROI_tChResp(rr,:,AnalysisBins));
            % get significant pvals
            Y = squeeze(data.pBinROI_tChResp(rr,:,AnalysisBins));
            Yt=reshape(mafdr(Y(:)),[nBands nBins]);
            X(Yt>info.pThr)=0;
            axes('position',[0.1 0.13 0.85 0.8])
           % h=contourfcmap(time,freqs,X,clev,y, 'lo',y2(1,:),'hi', y2(end,:),...
           %     'cbarloc','eastoutside','method','calccontour','evencb',1);
            h=contourfcmap(time,freqs,X,clev,y, 'lo',y2(1,:),'hi', y2(end,:),'method','calccontour');
            %              contourfcmap(time,freqs,X,clev,y, 'lo',[0.3 0.10 0.6],'hi', ...
            %                 [0.99 0.99 0.1],'method','calccontour')
            for i = 1:numel(h.l)
                h.l(i).Color='none';
            end
            
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            ytick = linspace(1.2,5.8,6);
            set(gca,'fontsize',25,'YTick',ytick,'ytickLabel',info.yticksLabels)
            set(gca,'XTick',info.xtick,'xtickLabel',info.xticklabel)
            info.fileName   = strcat(data.ROIs{rr},info.Analysis);
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
        end
        
        
    case 2
        % Plot Cluster activity
        figure(1); clf;
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,300]);
        info.Analysis   = [lock '_ClusterTVals' 'Tthr' num2str(opts.Kmeans.Thtr) 'K' num2str(opts.Kmeans.K)];

        load([dirPath '/Kmeans/' lock '_activity_Tthr' num2str(opts.Kmeans.Thtr) '_Kmeans' num2str(opts.Kmeans.K) ...
            'MultiBand' lock 'sublogPowernonLPCch']);

        freqs = 1:data.nBands; nBands = data.nBands;
        nBins = sum(data.AnalysisBins);
        time  = mean(data.Bins(data.AnalysisBins,:),2);
        
        for rr = 1:out.K
            X=squeeze(out.tBinCluster_tChResp(rr,:,data.AnalysisBins));
            Y=squeeze(out.pBinCluster_tChResp(rr,:,data.AnalysisBins));
            Yt=reshape(mafdr(Y(:)),[nBands nBins]);
            X(Yt>info.pThr)=0;            
            axes('position',[0.1 0.1 0.85 0.8])
            %contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            h=contourfcmap(time,freqs,X,clev,y, 'lo',y2(1,:),'hi', y2(end,:),'method','calccontour');
            %              contourfcmap(time,freqs,X,clev,y, 'lo',[0.3 0.10 0.6],'hi', ...
            %                 [0.99 0.99 0.1],'method','calccontour')
            for i = 1:numel(h.l)
                h.l(i).Color='none';
            end
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            set(gca,'fontsize',20,'ytickLabel',info.yticksLabels )
            set(gca,'XTick',info.xtick,'xtickLabel',info.xticklabel)
            info.fileName   = strcat('Cluster',num2str(rr),info.Analysis);
            %set(gca,'cLim',Zlims)
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
        end
        %fileName = [lock info.fileName];
        %renderChansSS2e(data.IDX,fileName)
        % semantic RTs correlations ROI
    case 3
        Zlims = [-3 3];
        clev = [Zlims(1):-1 1:Zlims(2)];
        y = brewermap(numel(clev)-1,'*RdBu');
        
        figure(1); clf;
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,250]);
        %colormap(cm);
        figure(2); clf;
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,400,250]);
        
        info.Analysis   = [lock 'StudyRTCorrROI_TVals'];
        freqs = 1:data.nBands; nBands = data.nBands;
        nBins = sum(data.AnalysisBins);
        time  = mean(data.Bins(data.AnalysisBins,:),2);
        time2 = [-0.25 0.25 0.75 1.25];
        xtick2 = [-0.25 0 0.25 0.5 0.75 1 1.25];
        xticklabel2 = {'','stim','','0.5','','1',''};
        for rr = 1:3
            chans=groups==rr;
            X = [];
            X2 = [];
            for ba=freqs
                [~,~,~,temp] = ttest(squeeze(data.dataToStudyRTsCorr(ba,chans,:)));
                X = [X;temp.tstat];
                [~,~,~,temp] = ttest(squeeze(data.StudyRTs_PrePostActModel_Ts3(ba,chans,:)));
                X2 = [X2;temp.tstat];
            end
                        
            figure(1);clf;
            axes('position',[0.13 0.15 0.85 0.8])
            %contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            h=contourfcmap(time,freqs,X,clev,y, 'lo',y2(1,:),'hi', y2(end,:),'method','calccontour');
            for i = 1:numel(h.l)
                h.l(i).Color='none';
            end
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            ytick = linspace(1.2,5.8,6);
            set(gca,'fontsize',30,'YTick',ytick,'ytickLabel',info.yticksLabels)  
            set(gca,'XTick',info.xtick,'xtickLabel',info.xticklabel)
            %set(gca,'cLim',Zlims)
            info.fileName   = strcat(data.ROIs{rr},info.Analysis);
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
            
            figure(2);clf;
            axes('position',[0.2 0.15 0.75 0.8])
            %axes('position',[0.1 0.13 0.85 0.8])
            %contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            h=contourfcmap(time2,freqs,X2,clev,y, 'lo',y2(1,:),'hi', y2(end,:),'method','calccontour');
            for i = 1:numel(h.l)
                h.l(i).Color='none';
            end
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            ytick = linspace(1.2,5.8,6);
            set(gca,'fontsize',30,'YTick',ytick,'ytickLabel',info.yticksLabels)  
            set(gca,'XTick',xtick2,'xtickLabel',xticklabel2)
            %set(gca,'cLim',Zlims)
            info.fileName   = strcat(data.ROIs{rr},'_bigbin_',info.Analysis);
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
            
        end
         
        %fileName = [lock 'StudyRTCorrClusterGroups'];
        %renderChansSS2e(data.KMeansStudy.IDX,fileName)
        
        % semantic RTs correlations Clusters
    case 4
        figure(1); clf;
        colormap(cm);
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,300]);
        
        info.Analysis   = [lock '_TVals'];
        freqs = 1:6;
        time  = mean(data.Bins,2);
        Zlims = [-3 3];
        for rr = 1:3
            X=squeeze(data.KMeansStudy.TTests(rr,:,:));
            X(~(abs(X)>info.thr))=0;
            info.fileName   = strcat('StudyRTCorrCluster',num2str(rr),info.Analysis);
            axes('position',[0.1 0.1 0.8 0.8])
            contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            set(gca,'fontsize',20,'ytickLabel',info.yticksLabels )
            set(gca,'cLim',Zlims)
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
        end
        fileName = [lock 'StudyRTCorrClusterGroups'];
        renderChansSS2e(data.KMeansStudy.IDX,fileName)
        
        % recogntion RTs Cluster Correlations
    case 5
        figure(1); clf;
        colormap(cm);
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,300]);
        
        info.Analysis   = [lock '_TVals'];
        freqs = 1:6;
        time  = mean(data.Bins,2);
        Zlims = [-3 3];
        for rr = 1:3
            X=squeeze(data.KMeansTest.TTests(rr,:,:));
            X(~(abs(X)>info.thr))=0;
            info.fileName   = strcat('TestRTCorrCluster',num2str(rr),info.Analysis);
            axes('position',[0.1 0.1 0.8 0.8])
            contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            set(gca,'fontsize',20,'ytickLabel',info.yticksLabels )
            set(gca,'cLim',Zlims)
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
        end
        fileName = [lock 'TestsRTCorrClusterGroups'];
        renderChansSS2e(data.KMeansTest.IDX,fileName)
    case 6
       Zlims = [-3 3];
        clev = [Zlims(1):-1 1:Zlims(2)];
        y = brewermap(numel(clev)-1,'*RdBu');
        
        figure(1); clf;
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,250]);
        figure(2); clf;
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,400,250]);
        %colormap(cm);
 
        info.Analysis   = [lock 'TestRTCorrROI_TVals'];
        freqs = 1:data.nBands; nBands = data.nBands;
        nBins = sum(data.AnalysisBins);
        time  = mean(data.Bins(data.AnalysisBins,:),2);
        time2 = [-0.25 0.25 0.75 1.25];
        xtick2 = [-0.25 0 0.25 0.5 0.75 1 1.25];
        xticklabel2 = {'','stim','','0.5','','1',''};
        for rr = 1:3
            chans=groups==rr;
            X = [];
            X2 = [];
            for ba=freqs
                [~,~,~,temp] = ttest(squeeze(data.dataToTestRTsCorr(ba,chans,:)));
                X = [X;temp.tstat];
                [~,~,~,temp] = ttest(squeeze(data.TestRTs_PrePostActModel_Ts3(ba,chans,:)));
                X2 = [X2;temp.tstat];
            end
            
            figure(1);clf;
            axes('position',[0.1 0.13 0.85 0.8])
            %contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            h=contourfcmap(time,freqs,X,clev,y, 'lo',y2(1,:),'hi', y2(end,:),'method','calccontour');
            for i = 1:numel(h.l)
                h.l(i).Color='none';
            end
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            ytick = linspace(1.2,5.8,6);
            set(gca,'fontsize',25,'YTick',ytick,'ytickLabel',info.yticksLabels)  
            set(gca,'XTick',info.xtick,'xtickLabel',info.xticklabel)
            %set(gca,'cLim',Zlims)
            info.fileName   = strcat(data.ROIs{rr},info.Analysis);
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
            
            figure(2);clf;
            axes('position',[0.15 0.13 0.8 0.8])
            %axes('position',[0.1 0.13 0.85 0.8])
            %contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            h=contourfcmap(time2,freqs,X2,clev,y, 'lo',y2(1,:),'hi', y2(end,:),'method','calccontour');
            for i = 1:numel(h.l)
                h.l(i).Color='none';
            end
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            ytick = linspace(1.2,5.8,6);
            set(gca,'fontsize',25,'YTick',ytick,'ytickLabel',info.yticksLabels)  
            set(gca,'XTick',xtick2,'xtickLabel',xticklabel2)
            %set(gca,'cLim',Zlims)
            info.fileName   = strcat(data.ROIs{rr},'_bigbin_',info.Analysis);
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
        end
end
end

function stats=univariateAnalysis(data,info)

% selection for statistics
X       = subSelectByCell(data.BinERP(subjs),info.cond1);
Y       = subSelectByCell(data.BinERP(subjs),info.cond2);
X       = concatenateCells(X);
Y       = concatenateCells(Y);
stats = BiVariateStatsWrapper(X,Y,info.groups);

info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
printStats(stats,info)

% selection for plotting
addpath Plotting/

X       = subSelectByCell(data.ERP(subjs),info.cond1);
Y       = subSelectByCell(data.ERP(subjs),info.cond2);

t = data.trialTime;
Z =[];
mean2 = @(X)(squeeze(nanmean(X,2)));
Z{1} = cell2mat(cellfun(mean2,X,'uniformoutput',0));
Z{2} = cell2mat(cellfun(mean2,Y,'uniformoutput',0));

info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
info.fileName   = strcat(info.Analysis,'AllChans');
plotWrapper(Z,t,info)

% by ROI
for iROI = 1:3
    chans = info.groups==iROI;
    Zroi =[];
    Zroi{1} = Z{1}(chans,:);
    Zroi{2} = Z{2}(chans,:);
    
    info.PVals      = stats.groupPVals(iROI,:);
    info.fileName   = sprintf('%s%s',info.Analysis,info.groupNames{iROI});
    plotWrapper(Zroi,t,info)
end

end

function stats=correlationAnalysis(data)


end