% Analysis Code for SS2 encoding data
%
% Analyses codes:
function SS2e_MBAnalysis(lock,AnalysisNum,opts)
inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';


addpath Analysis/
addpath Plotting/
addpath lib/

%dirPath         = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
dirPath         = ['~/Google Drive/Research/ECoG_SS2e/data_results/'] ;
%fileName        = ['allERSPs' band 'Group' lock 'sublogPowernonLPCCh'];
%fileName        = ['allMBAnalysis' lock 'sublogPowernLPClowEvokedVar'];
fileName        = ['allMBAnalysis' lock 'sublogPowernonLPCch'];

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
info.savePath   = ['~/Google Drive/Research/ECoG_SS2e' ...
    '/Plots/SS2e_MB_Analyses/a' num2str(AnalysisNum) '/'];
if ~exist(info.savePath,'dir')
    mkdir(info.savePath);
end

info.groupNames = {'IPS','SPL','AG'};

info.Bins       = data.Bins;
%info.yLimits    = [-1 1.6];
%info.yLimits    = [-4 2];
if strcmp(lock,'stim')
    info.xtick      = [0 0.4 0.8 1.2];
    info.toi        = [0 1.5];
    info.alpha      = 0.05/15;
    info.xticklabel = {'stim','0.4','0.8','1.2'};
elseif strcmp(lock,'RT')
    info.yAxisRightLoc = 1;
    info.xtick      = [-0.8 -0.4 0 0.4];
    info.xticklabel = {'-0.8','-0.4','resp','0.4'};
elseif strcmp(lock, 'preStim2')
    info.xtick      = [-.8 -0.4 0 0.4 0.8 1.2];
    info.toi        = [-1 1.5];
    info.alpha      = 0.05/26;
    info.xticklabel = {'-0.8','-0.4','stim','0.4','0.8','1.2'};
end

info.cols            = 'oc';
info.yticksLabels   = {'\delta','\theta','\alpha','\beta','l-\gamma','h-\gamma'};
info.thr            = 1;
info.pThr           = 0.05; % corrected.

% Analysis 1
clev = [-4:-1 1:4];
y = brewermap(numel(clev)-1,'*RdBu');
y2 = brewermap(1000,'*RdBu');
% x=parula(1200);
% d=downsample(1:size(x,1),floor(size(x,1)/(numel(clev)-2)));
% x=x(d,:);
% y=[x(1:numel(d)/2,:); 1 1 1; x(numel(d)/2+1:numel(d),:)];
%cm = [x(1:32,:);ones(32,3);x(33:64,:)];

%colormap(cm);
Zlims = [-4 4];
switch AnalysisNum
    % Plot ROIs level activity
    case 1
        figure(1); clf;
        %        colormap(cm);
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,300]);
        info.Analysis   = [lock '_SpectoROIsTVals'];
        
        freqs = 1:data.nBands; nBands = data.nBands;
        nBins = sum(data.AnalysisBins);
        time  = mean(data.Bins(data.AnalysisBins,:),2);
        for rr =1:3
            % select ROI spectrogram
            X = squeeze(data.tBinROI_tChResp(rr,:,data.AnalysisBins));
            % get significant pvals
            Y = squeeze(data.pBinROI_tChResp(rr,:,data.AnalysisBins));
            Yt=reshape(mafdr(Y(:)),[nBands nBins]);
            X(Yt>info.pThr)=0;
            axes('position',[0.1 0.1 0.85 0.8])
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
            set(gca,'fontsize',20,'YTick',freqs,'ytickLabel',info.yticksLabels)
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
        figure(1); clf;
        colormap(cm);
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,300]);
        
        info.Analysis   = [lock 'StudyRTCorrROI_TVals'];
        freqs = 1:6;
        time  = mean(data.Bins,2);
        for rr = 1:3
            chans=groups==rr;
            X = [];
            for ba=freqs
                [~,~,~,temp] = ttest(squeeze(data.dataToStudyRTsCorr(ba,chans,:)));
                X = [X;temp.tstat];
            end
            
            info.fileName   = strcat('StudyRTCorrROI',num2str(rr),info.Analysis);
            axes('position',[0.1 0.1 0.8 0.8])
            contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            set(gca,'fontsize',20,'ytickLabel',info.yticksLabels )
            set(gca,'cLim',Zlims)
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
        % recogntion RTs correlations ROIs
        figure(1); clf;
        colormap(cm);
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,300]);
        
        info.Analysis   = ['_TVals'];
        freqs = 1:6;
        time  = mean(data.Bins,2);
        Zlims = [-3 3];
        for rr = 1:3
            chans = data.ROIid==rr;
            X = [];
            for ba=freqs
                [~,~,~,temp] = ttest(squeeze(data.dataToTestRTsCorr(ba,chans,:)));
                X = [X;temp.tstat];
            end
            %X=squeeze(data.KMeansTest.TTests(rr,:,:));
            X(~(abs(X)>info.thr))=0;
            info.fileName   = strcat('TestRTCorrROIs',num2str(rr),info.Analysis);
            axes('position',[0.1 0.1 0.8 0.8])
            contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            set(gca,'fontsize',20,'ytickLabel',info.yticksLabels )
            set(gca,'cLim',Zlims)
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
        end
        %         %fileName = 'TestsRTCorrClusterGroups';
        %         %renderChansSS2e(data.KMeansTest.IDX,fileName)
        %     case 6
        %         % lasso recognition RT
        % %         for ss = 1:7
        % %             figure(1); clf;
        % %             set(gcf,'paperpositionmode','auto','color','white')
        % %             set(gcf,'position',[100,100,300,300]);
        % %             s=scatter(data.LassoTestRTs.Y{ss},data.LassoTestRTs.Yhat{ss},'ok');
        % %             set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',50)
        % %             l=lsline;
        % %             set(gca,'fontsize',20,'linewidth',4)
        % %             set(l,'linewidth',4)
        % %             xlabel(' RTs '); ylabel(' Prediction ')
        % %
        % %             xlim([0 3])
        % %             fN = [lock '_subj' num2str(ss) 'lassoRecogRTpred'];
        % %             fP = '~/Google Drive/Research/ECoG_SS2e/';
        % %             cPath = pwd;
        % %             cd(fP)
        % %             addpath(cPath)
        % %             addpath([cPath '/Plotting/'])
        % %
        % %             print(gcf,'-dsvg',fN)
        % %             eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
        % %             eval(['!rm ' fN '.svg'])
        % %         end
        %
        %         figure(2); clf; set(gcf, 'position', [100 100 500 400],'paperpositionmode','auto');
        %         sD = 200;
        %         shapes = 'ods><x+';
        %         % RTs (retrieval)
        %         axes('position',[0.2 0.1 0.40 0.85]); hold on;
        %         xlim([0 1]); ylim([-0.3 0.6])
        %         X = data.LassoTestRTs.corr;
        %         for ss= 1:7
        %             s=scatter(0.5+randn*0.1,X(ss),shapes(ss));
        %             set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
        %         end
        %         plot([0.2 0.8],[1 1]*mean(X),'linewidth',3,'color',0.4*ones(1,3))
        %         set(gca,'fontsize',20,'xtick',[],'ytick',[-0.5 0 0.5],'linewidth',2);
        %         ylabel(' Corr (r) ')
        %
        %         fN = [lock 'lassoRecogRTpred'];
        %         fP = '~/Google Drive/Research/ECoG_SS2e/';
        %         cPath = pwd;
        %         cd(fP)
        %         addpath(cPath)
        %         addpath([cPath '/Plotting/'])
        %
        %         print(gcf,'-dsvg',fN)
        %         eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
        %         eval(['!rm ' fN '.svg'])
        %
        %         cd(cPath)
        %     case 7
        %         for ss = 1:7
        %             figure(1); clf;
        %             set(gcf,'paperpositionmode','auto','color','white')
        %             set(gcf,'position',[100,100,300,300]);
        %             s=scatter(data.LassoStudyRTs.Y{ss},data.LassoStudyRTs.Yhat{ss},'ok');
        %             set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',50)
        %             l=lsline;
        %             set(gca,'fontsize',20,'linewidth',4)
        %             set(l,'linewidth',4)
        %             xlabel(' RTs '); ylabel(' Prediction ')
        %
        %             xlim([0 3])
        %             fN = [lock '_subj' num2str(ss) 'lassoStudyRTpred'];
        %             fP = '~/Google Drive/Research/ECoG_SS2e/';
        %             cPath = pwd;
        %             cd(fP)
        %             addpath(cPath)
        %             addpath([cPath '/Plotting/'])
        %
        %             print(gcf,'-dsvg',fN)
        %             eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
        %             eval(['!rm ' fN '.svg'])
        %         end
        %         figure(2); clf; set(gcf, 'position', [100 100 500 400],'paperpositionmode','auto');
        %         sD = 200;
        %         shapes = 'ods><x+';
        %         % RTs (retrieval)
        %         axes('position',[0.2 0.1 0.40 0.85]); hold on;
        %         xlim([0 1]); ylim([-0.3 0.6])
        %         X = data.LassoStudyRTs.corr;
        %         for ss= 1:7
        %             s=scatter(0.5+randn*0.1,X(ss),shapes(ss));
        %             set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
        %         end
        %         plot([0.2 0.8],[1 1]*mean(X),'linewidth',3,'color',0.4*ones(1,3))
        %         set(gca,'fontsize',20,'xtick',[],'ytick',[-0.5 0 0.5],'linewidth',2);
        %         ylabel(' Corr (r) ')
        %
        %         fN = [lock 'lassoStudyRTpred'];
        %         fP = '~/Google Drive/Research/ECoG_SS2e/';
        %         cPath = pwd;
        %         cd(fP)
        %         addpath(cPath)
        %         addpath([cPath '/Plotting/'])
        %
        %         print(gcf,'-dsvg',fN)
        %         eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
        %         eval(['!rm ' fN '.svg'])
        %
        %         cd(cPath)
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