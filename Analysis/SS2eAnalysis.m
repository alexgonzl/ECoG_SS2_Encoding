% Analysis Code for SS2 encoding data
%
% Analyses codes:
% 1) abstract vs concrete
% 2) subjective abs vs concrete
% 3) correct abs vs concrete
% 4) correct abs vs concrete, & correctly remembered
% 5) activity to test RT correlationa analysis
% 6) correlation using only correct abstract/concrete
% 7) correlation using only correct abstract/concrete & remembered
% 8) test RT correlation; broken down by encoding type (abs/conc)
% 9) test RT correlation; broken down by encoding type (remembered)
% 10) parametric slow-fast test RTs splits
% 11) encoding RT correlation; broken down by encoding type (remembered)
% 12) overall correlation to encoding RT
% 13) HGP time-courses without category split.        
% 14) enconding retrieval activity correlation hypothesis, higher SPL activity -> higher IPS activity for stronger memories

function stats=SS2eAnalysis(band,lock,AnalysisNum,hem)


addpath Analysis/
addpath Plotting/
addpath lib/

dirPath         = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
%fileName        = ['allERSPs' band 'Group' lock 'sublogPowernonLPCCh'];
fileName        = ['allERSPs' band 'Group' lock 'sublogPowernLPClowEvokedVar'];

dataPath = [dirPath 'group/Spectral_Data/'];
load([dataPath fileName '.mat'])

if ~strcmp(hem,'all')
    subjs = find(strcmp(data.subjHemsIds,hem));
else
    subjs = 1:numel(data.options.subjects);
end
nSubjs = numel(subjs);
subjChans = ismember(data.subjChans,subjs);

% channel groupings/ROIs
nChans = numel(data.ROIid);
groups = zeros(nChans,1);
groups(data.ROIid==1) = 1; %IPS
groups(data.ROIid==2) = 2; %SPL
groups(data.ROIid==3) = 3; %AG
groups(~subjChans)   = [];

correctAbsConc      = cell(nSubjs,1);
for ss = 1:nSubjs
    correctAbsConc{ss} = (data.conds{subjs(ss),1} & data.conds{subjs(ss),3}) |...
     (data.conds{subjs(ss),2} & data.conds{subjs(ss),4});
end
correctAbsConcRem   = cellfun(@and,correctAbsConc,data.conds(subjs,5),'UniformOutput',0);
%%
close all

info            =[];
info.groups     = groups;
info.rownames   = cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
info.groupNames = {'IPS','SPL','AG'};

info.Bins       = data.Bins;
info.alpha      = 0.05/15;
info.toi        = [0 1.5];
info.yLimits    = [-1 1.6];
%info.yLimits    = [-4 2];
if strcmp(lock,'stim')
    info.xtick      = [0 0.4 0.8 1.2];
    info.xticklabel = {'stim','0.4','0.8','1.2'};
elseif strcmp(lock,'RT')
    info.yAxisRightLoc = 1;
    info.xtick      = [-0.8 -0.4 0 0.4];
    info.xticklabel = {'-0.8','-0.4','resp','0.4'};
end
info.cols       = 'oc';

% Analysis 1
switch AnalysisNum
    case 1
        % Absstract vs Concrete
        info.Analysis   = [hem '_' lock band '_AbsConcByROI'];
        info.cond1      = data.conds(subjs,1); % abstract
        info.cond2      = data.conds(subjs,2); % concrete
        %info.legend 	= {'Abs','Conc'};
        
        stats=univariateAnalysis(data,info);
    case 2
        % Subjective Absstract vs Concrete
        info.Analysis   = [hem '_' lock band '_RespAbsConcByROI'];
        info.cond1           = data.conds(subjs,3); % abstract
        info.cond2           = data.conds(subjs,4); % concrete
        %info.legend 	= {'Resp. Abs','Resp. Conc'};
        
        stats=univariateAnalysis(data,info);
    case 3
        % Correct Absstract vs Concrete
        info.Analysis   = [hem '_' lock band '_CorrectAbsConcByROI'];
        info.cond1 	= cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
        info.cond2 	= cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
        %info.legend 	= {'Correct Abs','Correct Conc'};
        
        stats=univariateAnalysis(data,info);
    case 4
        % Correct Absstract vs Concrete & correctly rememberd
        % Correct Absstract vs Concrete
        info.Analysis   = [hem '_' lock band '_CorrectAbsConcRemByROI'];
        info.cond1 	= cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
        info.cond1 	= cellfun(@and,info.cond1,data.conds(subjs,5),'uniformoutput',0); % correctly remembered
        info.cond2 	= cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
        info.cond2 	= cellfun(@and,info.cond2,data.conds(subjs,5),'uniformoutput',0); % correct abstract
        %info.legend 	= {'Correct Abs','Correct Conc'};
        
        stats=univariateAnalysis(data,info);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % activity set 2: correlation to test RTs activity
    case 5
        % overall activity correlation to test RTs
        %info            =[];
        info.Analysis   = [hem '_' lock band '_testRTcorr'];
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend 	= info.groupNames;
        
        X = CorrByCell(data.BinERP(subjs),cellfun(@log10,data.testRTs(subjs),'uniformoutput',0));
        X = cellfun(@atanh,X,'uniformoutput',0); % convet to continous
        X = [X{:}]';
        
        stats = statsWrapperMat(X,groups);
        info.rownames 	= cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        %printStats(stats,info)
        
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.yLimits    = [-0.2 0.2];
        t = mean(data.Bins,2);
        cols = 'rbg';
        
        % correlation by region.
        for iROI = 1:3
            info.fileName 	= strcat(info.Analysis,info.groupNames{iROI});
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
        info.cols 		= 'rbg';
        info.yLimits    = [-5 3];
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.fileName 	= strcat(info.Analysis,'AllROIsScores');
        plotWrapper(Z,t,info)
    case 6
        % correlation using only correct abstract/concrete
        %info            =[];
        info.Analysis   = [hem '_' lock band '_correctAbsConc_testRTcorr'];
        
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend 	= info.groupNames;
        
        cond1 	= cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
        cond2 	= cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
        cond 	= cellfun(@or,cond1,cond2,'uniformoutput',0); % corrects
        
        X       = subSelectByCell(data.BinERP(subjs),cond);
        RTs     = cell(4,1);
        for ii = 1:4
            RTs{ii} = log10(data.testRTs{ii}(cond{ii}));
        end
        Y = CorrByCell(X,RTs);
        Y = cellfun(@atanh,Y,'uniformoutput',0); % convet to continous
        Y = [Y{:}]';
        
        stats = statsWrapperMat(Y,groups);
        info.rownames 	= cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        printStats(stats,info)
        
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.yLimits    = [-0.2 0.2];
        t = mean(data.Bins,2);
        cols = 'rbg';
        
        % correlation by region.
        for iROI = 1:3
            info.fileName 	= strcat(info.Analysis,info.groupNames{iROI});
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
        
        info.cols 		= 'rbg';
        info = rmfield(info,'PVals');
        info.cols 		= 'rbg';
        info.yLimits    = [-5 3];
        
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.fileName 	= strcat(info.Analysis,'AllROIsScores');
        plotWrapper(Z,t,info)        
    case 7
        % correlation using only correct abstract/concrete & remembered
        %info            =[];
        %info.Analysis   = [hem '_' lock band '_correctAbsConcRem_testRTcorr'];
        info.Analysis   = [hem '_' lock band '_validEnctrialsRem_testRTcorr'];
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend    = info.groupNames;
        
%         cond1   = cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
%         cond2   = cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
%         cond    = cellfun(@or,cond1,cond2,'uniformoutput',0);   % corrects
%         cond    = cellfun(@and,cond,data.conds(subjs,5),'uniformoutput',0);           % remembered
        cond    = data.conds(subjs,13);
        
        X       = subSelectByCell(data.BinERP(subjs),cond);
        RTs     = cell(nSubjs,1);
        for ss = 1:nSubjs
            RTs{ss} = log10(data.testRTs{subjs(ss)}(cond{ss}));
        end
        Y = CorrByCell(X,RTs);
        Y = cellfun(@atanh,Y,'uniformoutput',0); % convet to continous
        Y = [Y{:}]';
        
        stats = statsWrapperMat(Y,groups);
        info.rownames   = cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
%        printStats(stats,info)
        
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
        info.Analysis   = [hem '_' lock band '_correctAbsConc_testRTcorrByCond'];
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend 	= info.groupNames;
        
        cond1 	= cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
        cond2 	= cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
        
        X       = subSelectByCell(data.BinERP(subjs),cond1);
        Y       = subSelectByCell(data.BinERP(subjs),cond2);
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
        %info.rownames 	= cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        %info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        %printStats(stats,info)
        
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
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
        info.Analysis   = [hem '_' lock band '_correctAbsConcRem_testRTcorrByCond'];
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend    = info.groupNames;
        
        cond1   = cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
        cond1   = cellfun(@and,cond1,data.conds(subjs,5),'uniformoutput',0);
        cond2   = cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
        cond2   = cellfun(@and,cond2,data.conds(subjs,5),'uniformoutput',0);
        
        X       = subSelectByCell(data.BinERP(subjs),cond1);
        Y       = subSelectByCell(data.BinERP(subjs),cond2);
        
        RTsCond1     = cell(nSubjs,1);
        RTsCond2     = cell(nSubjs,1);
        for ss = 1:nSubjs
            RTsCond1{ss} = log10(data.testRTs{subjs(ss)}(cond1{ss}));
            RTsCond2{ss} = log10(data.testRTs{subjs(ss)}(cond2{ss}));
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
        info.yLimits    = [-0.25 0.4];
        
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
        info.Analysis   = [hem '_' lock band '_testRTSplitByQuantile'];
        info.groupNames = {'IPS','SPL','AG'};
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.legend     = 1;
        
        mean2 = @(X)(squeeze(nanmean(X,2)));
        nSplits = 4; nQuants = numel(data.RTquantiles);
        t = data.trialTime;
        Z =[];
        %         % plot all splits by region
        %         for iROI = 1:3
        %             chans = info.groups==iROI;
        %             for jj = 1:nSplits
        %                 cond1 = cell(nSubjs,nSplits);
        %                 cond2 = cell(nSubjs,nSplits);
        %                 for ss = 1:4
        %                     cond1{ss,jj} = data.testRTs{ss}<=data.testRTsQuants{ss}(jj);            % fast
        %                     cond2{ss,jj} = data.testRTs{ss}>=data.testRTsQuants{ss}(nQuants-jj+1);  % slow
        %                 end
        %
        %                 cond1(:,jj) = cellfun(@and,cond1(:,jj),correctAbsConcRem,'UniformOutput',0);
        %                 cond2(:,jj) = cellfun(@and,cond2(:,jj),correctAbsConcRem,'UniformOutput',0);
        %                 X       = subSelectByCell(data.ERP(subjs),cond1(:,jj));
        %                 Y       = subSelectByCell(data.ERP(subjs),cond2(:,jj));
        %                 temp    = cell2mat(cellfun(mean2,X,'uniformoutput',0));
        %                 Z{1}(jj,:) = nanmean(temp(chans,:));
        %                 temp = cell2mat(cellfun(mean2,Y,'uniformoutput',0));
        %                 Z{2}(jj,:) = nanmean(temp(chans,:));
        %             end
        %             info.fileName   = sprintf('%s%s',info.Analysis,info.groupNames{iROI});
        %             parametricTimeCoursePlot(Z,t,info)
        %         end
        
        %plot splits using hard boundaries
        RTBinsIds   = 1:2:nQuants;
        nRTBins     = numel(RTBinsIds);
        cond1 = cell(nSubjs,numel(RTBinsIds));
        for jj = 1:nRTBins
            for ss = 1:nSubjs
                if jj<=1
                    cond1{ss,jj} = data.testRTs{subjs(ss)}<=data.testRTsQuants{subjs(ss)}(1); % fast
                elseif jj>=nRTBins
                    cond1{ss,jj} = data.testRTs{subjs(ss)}>data.testRTsQuants{subjs(ss)}(nQuants); % slow
                else
                    cond1{ss,jj} = data.testRTs{subjs(ss)}<=data.testRTsQuants{subjs(ss)}(RTBinsIds(jj)+2); % fast
                    cond1{ss,jj} = (cond1{ss,jj}) & (data.testRTs{subjs(ss)}>data.testRTsQuants{subjs(ss)}(RTBinsIds(jj))); % fast
                end
            end
            cond1(:,jj) = cellfun(@and,cond1(:,jj),correctAbsConcRem,'uniformoutput',0);
        end
        
        for iROI = 1:3
            chans = info.groups==iROI;
            % plot all time course splits
            Z =[];
            for jj = 1:nRTBins
                X       = subSelectByCell(data.ERP(subjs),cond1(:,jj));
                temp    = cell2mat(cellfun(mean2,X,'uniformoutput',0));
                Z{1}(jj,:) = nanmean(temp(chans,:),1);
            end
            info.fileName       = sprintf('%s%s',info.Analysis,info.groupNames{iROI},'hardSplit');
            info.colID          = iROI;
            info.shade_toi      = 1;
            info.toi            = [0.5 0.7]; % time of interest in seconds for bar plot
            parametricTimeCoursePlot(Z,t,info)
            
            % bar plots for time of interest
            Z =[];
            info2               = info;
            samps = data.trialTime>=info2.toi(1) & data.trialTime<=info2.toi(2);
            info2.yLimits       = [-1 1.6];
            info2.yTicks         = [-.5 0 0.75 1.5];
            info2.aspectRatio   = [400 300];
            info2.xTicks        = [1 nRTBins];
            info2.XTickLabel    = {'Fastest','Slowest'};
            if iROI==1
                info2.colors        =  [linspace(0.5,1.0,5)',linspace(0.0,0.8,5)',linspace(0.0,0.8,5)'];
            elseif iROI==2
                info2.colors        =  [linspace(0,0.8,5)',linspace(0.0,0.8,5)',linspace(0.5,1.0,5)'];
            else
                info2.colors        =  [linspace(0,0.8,5)',linspace(0.5,1.0,5)',linspace(0.0,0.8,5)'];
            end
            info2.save          = 1;
            for jj = 1:nRTBins
                X       = subSelectByCell(data.ERP(subjs),cond1(:,jj));
                temp    = cell2mat(cellfun(mean2,X,'uniformoutput',0));
                Z{jj} = nanmean(temp(chans,samps),2);
            end
            toiStr           = strrep(strcat(num2str(info2.toi(1)),'-to-',num2str(info2.toi(2))),'.','p');
            info2.fileName   = strcat(info.Analysis,info.groupNames{iROI},'barQuintiles',toiStr);
            barPlotWithErrors(Z,info2);
        end
        
        %         Z =[];
        %         % plot by difference in quantile
        %         for iROI = 1:3
        %             chans = info.groups==iROI;
        %             for jj = 1:nSplits
        %                 cond1 = cell(nSubjs,nSplits);
        %                 cond2 = cell(nSubjs,nSplits);
        %                 for ss = 1:4
        %                     cond1{ss,jj} = data.testRTs{ss}<=data.testRTsQuants{ss}(jj); % fast
        %                     cond2{ss,jj} = data.testRTs{ss}>=data.testRTsQuants{ss}(nQuants-jj+1); % slow
        %                 end
        %                 X       = subSelectByCell(data.ERP(subjs),cond1(:,jj));
        %                 Y       = subSelectByCell(data.ERP(subjs),cond2(:,jj));
        %                 temp1    = cell2mat(cellfun(mean2,X,'uniformoutput',0));
        %                 temp2 = cell2mat(cellfun(mean2,Y,'uniformoutput',0));
        %                 Z{1}(jj,:) = nanmean(temp1(chans,:))-nanmean(temp2(chans,:));
        %             end
        %             info.fileName   = sprintf('%s%s',info.Analysis,info.groupNames{iROI},'SplitDiff');
        %             parametricTimeCoursePlot(Z,t,info);
        %         end
        %
        %         info = rmfield(info,'legend');
        %         % plot by split
        %         for jj = 1:nSplits
        %             cond1 = cell(nSubjs,1);
        %             cond2 = cell(nSubjs,1);
        %             for ss = 1:4
        %                 cond1{ss} = data.testRTs{ss}<=data.testRTsQuants{ss}(jj); % fast
        %                 cond2{ss} = data.testRTs{ss}>=data.testRTsQuants{ss}(nQuants-jj+1); % slow
        %             end
        %             X       = subSelectByCell(data.ERP(subjs),cond1);
        %             Y       = subSelectByCell(data.ERP(subjs),cond2);
        %
        %             t = data.trialTime;
        %             Z =[];
        %             mean2 = @(X)(squeeze(nanmean(X,2)));
        %             Z{1} = cell2mat(cellfun(mean2,X,'uniformoutput',0));
        %             Z{2} = cell2mat(cellfun(mean2,Y,'uniformoutput',0));
        %
        %             info.cols = 'ly';
        %             info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        %             % by ROI
        %             for iROI = 1:3
        %                 chans = info.groups==iROI;
        %                 Zroi =[];
        %                 Zroi{1} = Z{1}(chans,:);
        %                 Zroi{2} = Z{2}(chans,:);
        %
        %                 info.fileName   = sprintf('%s%s',info.Analysis,info.groupNames{iROI},'split',num2str(jj));
        %                 plotWrapper(Zroi,t,info)
        %             end
        %         end
    case 11
        % Analyze activity - encoding RT correlation; broken down by encoding type (remembered)
        info.Analysis   = [hem '_' lock band '_correctAbsConcRem_studyRTcorrByCond'];
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend    = info.groupNames;
        
        cond1   = cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
        cond1   = cellfun(@and,cond1,data.conds(subjs,5),'uniformoutput',0);
        cond2   = cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
        cond2   = cellfun(@and,cond2,data.conds(subjs,5),'uniformoutput',0);
        
        X       = subSelectByCell(data.BinERP(subjs),cond1);
        Y       = subSelectByCell(data.BinERP(subjs),cond2);
        RTsCond1     = cell(4,1);
        RTsCond2     = cell(4,1);
        for ii = 1:4
            RTsCond1{ii} = log10(data.studyRTs{ii}(cond1{ii}));
            RTsCond2{ii} = log10(data.studyRTs{ii}(cond2{ii}));
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
        info.yLimits    = [-0.3 0.45];
        
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
    case 12
        % overall correlation to encoding RT
        %info            =[];
        info.Analysis   = [hem '_' lock band '_correctAbsConcRem_encRTcorr'];
        info.groupNames = {'IPS','SPL','AG'};
        %info.legend    = info.groupNames;
        
        cond1   = cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
        cond2   = cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
        cond    = cellfun(@or,cond1,cond2,'uniformoutput',0);   % corrects
        cond    = cellfun(@and,cond,data.conds(subjs,5),'uniformoutput',0);           % remembered
        
        X       = subSelectByCell(data.BinERP(subjs),cond);
        RTs     = cell(nSubjs,1);
        for ii = 1:nSubjs
            RTs{ii} = log10(data.studyRTs{subjs(ii)}(cond{ii}));
        end
        Y = CorrByCell(X,RTs);
        Y = cellfun(@atanh,Y,'uniformoutput',0); % convet to continous
        Y = [Y{:}]';
        
        stats = statsWrapperMat(Y,groups);
        info.rownames   = cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        %info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        %printStats(stats,info)
        
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
        
    case 13
        %HGP time-courses without category split.        
        info.Analysis   = [hem '_' lock band '_CorrectAllByROI'];
        info.cond1 	= cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
        info.cond1 	= cellfun(@and,info.cond1,data.conds(subjs,5),'uniformoutput',0); % correctly remembered
        info.cond2 	= cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
        info.cond2 	= cellfun(@and,info.cond2,data.conds(subjs,5),'uniformoutput',0); % correct remembered
        info.cond3  = cellfun(@or, info.cond1,info.cond2,'uniformoutput',0); %
        cond        = data.conds(subjs,13);
        %info.legend 	= {'Correct Abs','Correct Conc'};
        
        
        X       = subSelectByCell(data.BinERP(subjs),cond);
        X       = concatenateCells(X);
        stats = statsWrapper(X,info.groups);
        %info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        %printStats(stats,info);
         
        X       = subSelectByCell(data.ERP(subjs),cond);
        X       = concatenateCells(X);        
        
        t = data.trialTime;
        Z =[];
        mean1 = @(X)(squeeze(nanmean(X,1)));
        Z{1} = cell2mat(cellfun(mean1,X,'uniformoutput',0));        
        info.savePath   = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        
        %info.alpha  = 0.05/15;
        cols        = 'rbg';
        % by ROI
        for iROI = 1:3
            chans = info.groups==iROI;
            Zroi =[];
            Zroi{1} = Z{1}(chans,:);            
            
            info.PVals      = stats.groupPVals(iROI,:);
            info.fileName   = sprintf('%s%s',info.Analysis,info.groupNames{iROI});
            info.cols       = cols(iROI);
            plotWrapper(Zroi,t,info)
        end
    case 14
        % enconding retrieval activity correlation
        % hypothesis, higher SPL activity -> higher IPS activity for stronger memories

        % load retrieval data.
        ret = load('/Users/alexandergonzalez/Documents/ECOG_Mem_Encoding/Data/allERSPshgamGroupstimLocksublogPowernonLPCleasL1TvalCh10.mat');
        % exclude subject 6, different expt.
        ret.data.BinERP(6) =[];
        Y = subSelectByCell(ret.data.BinERP(subjs),data.studyIDatTest);
        
        info.Analysis   = [hem '_' lock band '_enc_stimHGPret'];
        info.cond1  = cellfun(@and,data.conds(subjs,1),data.conds(subjs,3),'uniformoutput',0); % correct abstract
        info.cond1  = cellfun(@and,info.cond1,data.conds(subjs,5),'uniformoutput',0); % correctly remembered
        info.cond2  = cellfun(@and,data.conds(subjs,2),data.conds(subjs,4),'uniformoutput',0); % correct concrete
        info.cond2  = cellfun(@and,info.cond2,data.conds(subjs,5),'uniformoutput',0); % correct remembered
        info.cond3  = cellfun(@or, info.cond1,info.cond2,'uniformoutput',0); %

        X       = subSelectByCell(data.BinERP(subjs),info.cond3);
        X       = concatenateCells(X);
        Y       = subSelectByCell(Y,info.cond3);
        Y       = concatenateCells(Y);
        
        Z       = [];
        cnt     = 1;
        if strcmp(hem, 'l')
            for ss = [1 3 4]
                chans = find(data.subjChans==ss);                
                IPSchans    = find(groups(chans)==1)';
                SPLchans    = find(groups(chans)==2)';
                for ii = IPSchans
                    for jj = SPLchans
                        Z{cnt} = corr(X{jj}(:,7),Y{ii});
                        cnt = cnt +1;
                    end
                end
            end            
        else
            fprintf('not ready')
            return
        end
        %  channel correlation from encoding to retrieval
        Z = cellfun(@corr,X,Y,'uniformoutput',0);
        ZM = permute(reshape(cell2mat(Z),[17 nChans 12]),[2 1 3]);
        
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