
% Analysis Code for SS2 encoding data, stimulus locked, high gamma power.

function stats=rt_ThetaAnalysis(AnalysisNum)

addpath Analysis/
addpath lib/

dirPath       = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
fileName = 'allERSPsthetaGroupRTLocksublogPowernonLPCCh';

dataPath = [dirPath 'group/Spectral_Data/'];
load([dataPath fileName '.mat'])

% channel groupings/ROIs
nChans = numel(data.ROIid);
groups = zeros(nChans,1);
groups(data.ROIid==1) = 1; %IPS
groups(data.ROIid==2) = 2; %SPL
groups(data.ROIid==3) = 3; %AG


%%

info            =[];
info.groups     = groups;
info.rownames   = cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
info.groupNames = {'IPS','SPL','AG'};

info.Bins       = data.Bins;
info.alpha      = 0.01;
info.xtick      = [-0.8 -0.4 0 0.4];
info.yLimits    = [-6 2];
info.xticklabel = {'-0.8','-0.4','resp','0.4'};
info.cols       = 'oc';
info.yAxisRightLoc = 1;

switch AnalysisNum
    % Analysis 1.
    
    case 1
        % Absstract vs Concrete
        info.Analysis   = 'RT_THP_AbsConcByROI';
        info.cond1           = data.conds(:,1); % abstract
        info.cond2           = data.conds(:,2); % concrete
        %info.legend 	= {'Abs','Conc'};
        
        stats=univariateAnalysis(data,info);
        
    case 2
        % Subjective Absstract vs Concrete
        info.Analysis   = 'RT_THP_RespAbsConcByROI';
        info.cond1           = data.conds(:,3); % abstract
        info.cond2           = data.conds(:,4); % concrete
        %info.legend 	= {'Resp. Abs','Resp. Conc'};
        
        stats=univariateAnalysis(data,info);
        
    case 3
        % Correct Absstract vs Concrete
        info.Analysis   = 'RT_THP_CorrectAbsConcByROI';
        info.cond1 	= cellfun(@and,data.conds(:,1),data.conds(:,3),'uniformoutput',0); % correct abstract
        info.cond2 	= cellfun(@and,data.conds(:,2),data.conds(:,4),'uniformoutput',0); % correct concrete
        %info.legend 	= {'Correct Abs','Correct Conc'};
        
        stats=univariateAnalysis(data,info);
    case 4
        % Correct Absstract vs Concrete & correctly rememberd
        info.Analysis   = 'rtTHP_CorrectAbsConcRemByROI';
        info.cond1 	= cellfun(@and,data.conds(:,1),data.conds(:,3),'uniformoutput',0); % correct abstract
        info.cond1 	= cellfun(@and,info.cond1,data.conds(:,5),'uniformoutput',0); % correctly remembered
        info.cond2 	= cellfun(@and,data.conds(:,2),data.conds(:,4),'uniformoutput',0); % correct concrete
        info.cond2 	= cellfun(@and,info.cond2,data.conds(:,5),'uniformoutput',0); % correct abstract
        %info.legend 	= {'Correct Abs','Correct Conc'};
        
        stats=univariateAnalysis(data,info);

    case 5
        % overall correlation
        %info            =[];
        info.Analysis   = 'rtTHP_testRTcorr';
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend    = info.groupNames;
        
        X = CorrByCell(data.BinERP,cellfun(@log10,data.testRTs,'uniformoutput',0));
        X = cellfun(@atanh,X,'uniformoutput',0); % convet to continous
        X = [X{:}]';
        
        stats = statsWrapperMat(X,groups);
        info.rownames   = cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        printStats(stats,info)
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.yLimits    = [-0.2 0.2];
        t = mean(data.Bins,2);
        cols = 'rbg';
        
        % correlation by region.
        for iROI = 1:3
            info.fileName   = strcat(info.Analysis,info.groupNames{iROI});
            info.cols       = cols(iROI);
            info.PVals      = stats.groupPVals(iROI,:);
            Z =[];
            Z{1} = X(groups==iROI,:);
            plotWrapper(Z,t,info)
        end
        
        
        t = mean(data.Bins,2);
        Z =[];
        Z{1} = repmat(stats.groupScores(1,:),[3,1]); % hard coded for now.
        Z{2} = repmat(stats.groupScores(2,:),[3,1]);
        Z{3} = repmat(stats.groupScores(3,:),[3,1]);
        
        info = rmfield(info,'PVals');
        info.cols       = 'rbg';
        info.yLimits    = [-5 3];
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.fileName   = strcat(info.Analysis,'AllROIsScores');
        plotWrapper(Z,t,info)
    case 6
        %% correlation using only correct abstract/concrete
        %info            =[];
        info.Analysis   = 'rtTHPcorrectAbsConc_testRTcorr';
        
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend    = info.groupNames;
        
        cond1   = cellfun(@and,data.conds(:,1),data.conds(:,3),'uniformoutput',0); % correct abstract
        cond2   = cellfun(@and,data.conds(:,2),data.conds(:,4),'uniformoutput',0); % correct concrete
        cond    = cellfun(@or,cond1,cond2,'uniformoutput',0); % corrects
        
        X       = subSelectByCell(data.BinERP,cond);
        RTs     = cell(4,1);
        for ii = 1:4
            RTs{ii} = log10(data.testRTs{ii}(cond{ii}));
        end
        Y = CorrByCell(X,RTs);
        Y = cellfun(@atanh,Y,'uniformoutput',0); % convet to continous
        Y = [Y{:}]';
        
        stats = statsWrapperMat(Y,groups);
        info.rownames   = cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        printStats(stats,info)
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.yLimits    = [-0.2 0.2];
        t = mean(data.Bins,2);
        cols = 'rbg';
        
        % correlation by region.
        for iROI = 1:3
            info.fileName   = strcat(info.Analysis,info.groupNames{iROI});
            info.cols       = cols(iROI);
            info.PVals      = stats.groupPVals(iROI,:);
            Z =[];
            Z{1} = Y(groups==iROI,:);
            plotWrapper(Z,t,info)
        end
        
        t = mean(data.Bins,2);
        Z =[];
        Z{1} = repmat(stats.groupScores(1,:),[3,1]); % hard coded for now.
        Z{2} = repmat(stats.groupScores(2,:),[3,1]);
        Z{3} = repmat(stats.groupScores(3,:),[3,1]);
        
        info.cols       = 'rbg';
        info = rmfield(info,'PVals');
        info.cols       = 'rbg';
        info.yLimits    = [-5 3];
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.fileName   = strcat(info.Analysis,'AllROIsScores');
        plotWrapper(Z,t,info)
        
    case 7 
         % correlation using only correct abstract/concrete & remembered
        info.Analysis   = 'rtTHPcorrectAbsConcRem_testRTcorr';        
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend    = info.groupNames;
        
        cond1   = cellfun(@and,data.conds(:,1),data.conds(:,3),'uniformoutput',0); % correct abstract
        cond2   = cellfun(@and,data.conds(:,2),data.conds(:,4),'uniformoutput',0); % correct concrete
        cond    = cellfun(@or,cond1,cond2,'uniformoutput',0);   % corrects
        cond    = cellfun(@and,cond,data.conds(:,5),'uniformoutput',0);           % remembered
        
        X       = subSelectByCell(data.BinERP,cond);
        RTs     = cell(4,1);
        for ii = 1:4
            RTs{ii} = log10(data.testRTs{ii}(cond{ii}));
        end
        Y = CorrByCell(X,RTs);
        Y = cellfun(@atanh,Y,'uniformoutput',0); % convet to continous
        Y = [Y{:}]';
        
        stats = statsWrapperMat(Y,groups);
        info.rownames   = cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        printStats(stats,info)
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.yLimits    = [-0.2 0.2];
        t = mean(data.Bins,2);
        cols = 'rbg';
        
        % correlation by region.
        for iROI = 1:3
            info.fileName   = strcat(info.Analysis,info.groupNames{iROI});
            info.cols       = cols(iROI);
            info.PVals      = stats.groupPVals(iROI,:);
            Z =[];
            Z{1} = Y(groups==iROI,:);
            plotWrapper(Z,t,info)
        end
        
        t = mean(data.Bins,2);
        Z =[];
        Z{1} = repmat(stats.groupScores(1,:),[3,1]); % hard coded for now.
        Z{2} = repmat(stats.groupScores(2,:),[3,1]);
        Z{3} = repmat(stats.groupScores(3,:),[3,1]);
        
        info.cols       = 'rbg';
        info = rmfield(info,'PVals');
        info.cols       = 'rbg';
        info.yLimits    = [-5 3];
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.fileName   = strcat(info.Analysis,'AllROIsScores');
        plotWrapper(Z,t,info)
    case 8
        % Analyze activity - test RT correlation; broken down by encoding
        % type
        info.Analysis   = 'rtTHPcorrectAbsConc_testRTcorrByCond';
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend    = info.groupNames;
        
        cond1   = cellfun(@and,data.conds(:,1),data.conds(:,3),'uniformoutput',0); % correct abstract
        cond2   = cellfun(@and,data.conds(:,2),data.conds(:,4),'uniformoutput',0); % correct concrete
        
        X       = subSelectByCell(data.BinERP,cond1);
        Y       = subSelectByCell(data.BinERP,cond2);
        RTsCond1     = cell(4,1);
        RTsCond2     = cell(4,1);
        for ii = 1:4
            RTsCond1{ii} = log10(data.testRTs{ii}(cond1{ii}));
            RTsCond2{ii} = log10(data.testRTs{ii}(cond2{ii}));
        end
        Xc = CorrByCell(X,RTsCond1);
        Yc = CorrByCell(Y,RTsCond2);
        Xc = cellfun(@atanh,Xc,'uniformoutput',0); % convet to continous
        Yc = cellfun(@atanh,Yc,'uniformoutput',0); % convet to continous
        Xc = [Xc{:}]' ; Yc = [Yc{:}]';
        
        %stats = BiVariateStatsWrapper(Xc,Yc,info.groups);
        %info.rownames  = cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        %info.savePath  = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        %printStats(stats,info)
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.yLimits    = [-0.4 0.4];
        
        t = mean(data.Bins,2);
        Z =[];
        Z{1} = Xc;
        Z{2} = Yc;
        
        % by ROI
        for iROI = 1:3
            chans = info.groups==iROI;
            Zroi =[];
            Zroi{1} = Z{1}(chans,:);
            Zroi{2} = Z{2}(chans,:);
            
            [~,info.PVals]      = ttest(Zroi{1},Zroi{2});
            %info.PVals      = stats.groupPVals(iROI,:);
            info.fileName   = sprintf('%s%s',info.Analysis,info.groupNames{iROI});
            plotWrapper(Zroi,t,info)
        end
    case 9
         % Analyze activity - test RT correlation; broken down by encoding type (remembered)
        info.Analysis   = 'rtTHPcorrectAbsConcRem_testRTcorrByCond';
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend    = info.groupNames;
        
        cond1   = cellfun(@and,data.conds(:,1),data.conds(:,3),'uniformoutput',0); % correct abstract
        cond1   = cellfun(@and,cond1,data.conds(:,5),'uniformoutput',0);
        cond2   = cellfun(@and,data.conds(:,2),data.conds(:,4),'uniformoutput',0); % correct concrete
        cond2   = cellfun(@and,cond2,data.conds(:,5),'uniformoutput',0);

        X       = subSelectByCell(data.BinERP,cond1);
        Y       = subSelectByCell(data.BinERP,cond2);
        RTsCond1     = cell(4,1);
        RTsCond2     = cell(4,1);
        for ii = 1:4
            RTsCond1{ii} = log10(data.testRTs{ii}(cond1{ii}));
            RTsCond2{ii} = log10(data.testRTs{ii}(cond2{ii}));
        end
        Xc = CorrByCell(X,RTsCond1);
        Yc = CorrByCell(Y,RTsCond2);
        Xc = cellfun(@atanh,Xc,'uniformoutput',0); % convet to continous
        Yc = cellfun(@atanh,Yc,'uniformoutput',0); % convet to continous
        Xc = [Xc{:}]' ; Yc = [Yc{:}]';
        
        %stats = BiVariateStatsWrapper(Xc,Yc,info.groups);
        %info.rownames  = cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        %info.savePath  = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        %printStats(stats,info)
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.yLimits    = [-0.4 0.4];
        
        t = mean(data.Bins,2);
        Z =[];        
        Z{1} = Xc;
        Z{2} = Yc;
        
        % by ROI
        for iROI = 1:3
            chans = info.groups==iROI;
            Zroi =[];
            Zroi{1} = Z{1}(chans,:);
            Zroi{2} = Z{2}(chans,:);
            
            [~,info.PVals]      = ttest(Zroi{1},Zroi{2});
            %info.PVals      = stats.groupPVals(iROI,:);
            info.fileName   = sprintf('%s%s',info.Analysis,info.groupNames{iROI});
            plotWrapper(Zroi,t,info)
        end
    case 10
        % parametric slow-fast test RTs splits
end
end
function stats=univariateAnalysis(data,info)

% selection for statistics
X       = subSelectByCell(data.BinERP,info.cond1);
Y       = subSelectByCell(data.BinERP,info.cond2);
X       = concatenateCells(X);
Y       = concatenateCells(Y);
stats = BiVariateStatsWrapper(X,Y,info.groups);

info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
printStats(stats,info)

% selection for plotting
addpath Plotting/

X       = subSelectByCell(data.ERP,info.cond1);
Y       = subSelectByCell(data.ERP,info.cond2);

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
