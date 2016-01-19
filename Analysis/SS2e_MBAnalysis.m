% Analysis Code for SS2 encoding data
%
% Analyses codes:

function SS2e_MBAnalysis(lock,AnalysisNum)
inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';


addpath Analysis/
addpath Plotting/
addpath lib/

%dirPath         = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
dirPath         = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/data_results/';
%fileName        = ['allERSPs' band 'Group' lock 'sublogPowernonLPCCh'];
fileName        = ['allMBAnalysis' lock 'sublogPowernLPClowEvokedVar'];

load([dirPath fileName '.mat'])

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
info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
info.groupNames = {'IPS','SPL','AG'};

info.Bins       = data.Bins;
info.alpha      = 0.05/15;
info.toi        = [0 1.5];
%info.yLimits    = [-1 1.6];
%info.yLimits    = [-4 2];
if strcmp(lock,'stim')
    info.xtick      = [0 0.4 0.8 1.2];
    info.xticklabel = {'stim','0.4','0.8','1.2'};
elseif strcmp(lock,'RT')
    info.yAxisRightLoc = 1;
    info.xtick      = [-0.8 -0.4 0 0.4];
    info.xticklabel = {'-0.8','-0.4','resp','0.4'};
end
info.cols            = 'oc';
info.yticksLabels    = {'\delta','\theta','\alpha','\beta','l-\gamma','h-\gamma'};
info.thr            = 1;
% Analysis 1
x=colormap;
cm = [x(1:32,:);ones(32,3);x(33:64,:)];
colormap(cm);
switch AnalysisNum
    % Plot ROIs level activity
    case 1
        figure(1); clf;
        colormap(cm);
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,300]);
        
        
        info.Analysis   = [lock '_SpectoROIsTVals'];
        
        freqs = 1:6;
        time  = mean(data.Bins,2);
        Zlims = [-3 3];
        for rr =1:3
            X = squeeze(data.BinROITTests(rr,:,:));
            X(~(abs(X)>info.thr))=0;
            info.fileName   = strcat(data.ROIs{rr},info.Analysis);
            %X(~(abs(X)>info.thr))=0;
            axes('position',[0.1 0.1 0.8 0.8])
            contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            set(gca,'fontsize',20,'ytickLabel',info.yticksLabels )
            set(gca,'cLim',Zlims)
            
            %             colorbar
            %             %axes('position',[0.95 yPos(2) 0.05 height])
            %             cb=colorbar;
            %             set(cb,'position',[0.85 0.3 0.05 0.3])
            %             %set(cb,'CLim',Zlims)
            %             set(cb,'ytick',[Zlims(1) 0 Zlims(2)],'yticklabel',[Zlims(1) 0 Zlims(2)])
            %             %axis off
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
        end
        
    case 2
        % Plot Cluster activity
        figure(1); clf;
        colormap(cm);
        set(gcf,'paperpositionmode','auto','color','white')
        set(gcf,'position',[100,100,500,300]);
        
        info.Analysis   = [lock '_TVals'];
        freqs = 1:6;
        time  = mean(data.Bins,2);
        Zlims = [-3 3];
        for rr = 1:3
            X=squeeze(data.KMeansERA.TTests(rr,:,:));
            X(~(abs(X)>info.thr))=0;
            info.fileName   = strcat('Cluster',num2str(rr),info.Analysis);
            axes('position',[0.1 0.1 0.8 0.8])
            contourf(time,freqs,X,[-10 -3 -2 0 2 3 10]);
            hold on;
            plot([0 0],[1 6],'k','linewidth',4)
            set(gca,'fontsize',20,'ytickLabel',info.yticksLabels )
            set(gca,'cLim',Zlims)
            print(gcf, '-dtiff','-r300', [info.savePath info.fileName])
        end
        fileName = [lock 'ERAClusterGroups'];
        renderChansSS2e(data.KMeansERA.IDX,fileName)
    case 3
        % semantic RTs correlations
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
    case 4
        % recogntion RTs Cluster Correlations
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
    case 5
        % recogntion RTs correlations
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
        %fileName = 'TestsRTCorrClusterGroups';
        %renderChansSS2e(data.KMeansTest.IDX,fileName)
    case 6
        % lasso recognition RT
%         for ss = 1:7
%             figure(1); clf;
%             set(gcf,'paperpositionmode','auto','color','white')
%             set(gcf,'position',[100,100,300,300]);
%             s=scatter(data.LassoTestRTs.Y{ss},data.LassoTestRTs.Yhat{ss},'ok');
%             set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',50)
%             l=lsline;
%             set(gca,'fontsize',20,'linewidth',4)
%             set(l,'linewidth',4)
%             xlabel(' RTs '); ylabel(' Prediction ')
%             
%             xlim([0 3])
%             fN = [lock '_subj' num2str(ss) 'lassoRecogRTpred'];
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
        
        figure(2); clf; set(gcf, 'position', [100 100 500 400],'paperpositionmode','auto');
        sD = 200;
        shapes = 'ods><x+';
        % RTs (retrieval)
        axes('position',[0.2 0.1 0.40 0.85]); hold on;
        xlim([0 1]); ylim([-0.3 0.6])
        X = data.LassoTestRTs.corr;
        for ss= 1:7
            s=scatter(0.5+randn*0.1,X(ss),shapes(ss));
            set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
        end
        plot([0.2 0.8],[1 1]*mean(X),'linewidth',3,'color',0.4*ones(1,3))
        set(gca,'fontsize',20,'xtick',[],'ytick',[-0.5 0 0.5],'linewidth',2);
        ylabel(' Corr (r) ')
        
        fN = [lock 'lassoRecogRTpred'];
        fP = '~/Google Drive/Research/ECoG_SS2e/';
        cPath = pwd;
        cd(fP)
        addpath(cPath)
        addpath([cPath '/Plotting/'])
        
        print(gcf,'-dsvg',fN)
        eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
        eval(['!rm ' fN '.svg'])
        
        cd(cPath)
    case 7
        for ss = 1:7
            figure(1); clf;
            set(gcf,'paperpositionmode','auto','color','white')
            set(gcf,'position',[100,100,300,300]);
            s=scatter(data.LassoStudyRTs.Y{ss},data.LassoStudyRTs.Yhat{ss},'ok');
            set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',50)
            l=lsline;
            set(gca,'fontsize',20,'linewidth',4)
            set(l,'linewidth',4)
            xlabel(' RTs '); ylabel(' Prediction ')
            
            xlim([0 3])
            fN = [lock '_subj' num2str(ss) 'lassoStudyRTpred'];
            fP = '~/Google Drive/Research/ECoG_SS2e/';
            cPath = pwd;
            cd(fP)
            addpath(cPath)
            addpath([cPath '/Plotting/'])
            
            print(gcf,'-dsvg',fN)
            eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
            eval(['!rm ' fN '.svg'])
        end
        figure(2); clf; set(gcf, 'position', [100 100 500 400],'paperpositionmode','auto');
        sD = 200;
        shapes = 'ods><x+';
        % RTs (retrieval)
        axes('position',[0.2 0.1 0.40 0.85]); hold on;
        xlim([0 1]); ylim([-0.3 0.6])
        X = data.LassoStudyRTs.corr;
        for ss= 1:7
            s=scatter(0.5+randn*0.1,X(ss),shapes(ss));
            set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD,'lineWidth',4)
        end
        plot([0.2 0.8],[1 1]*mean(X),'linewidth',3,'color',0.4*ones(1,3))
        set(gca,'fontsize',20,'xtick',[],'ytick',[-0.5 0 0.5],'linewidth',2);
        ylabel(' Corr (r) ')
        
        fN = [lock 'lassoStudyRTpred'];
        fP = '~/Google Drive/Research/ECoG_SS2e/';
        cPath = pwd;
        cd(fP)
        addpath(cPath)
        addpath([cPath '/Plotting/'])
        
        print(gcf,'-dsvg',fN)
        eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
        eval(['!rm ' fN '.svg'])
        
        cd(cPath)
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