% Analysis Code for SS2 encoding data, stimulus locked, high gamma power.


addpath Analysis/
addpath lib/

dirPath       = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
fileName = 'allERSPshgamGroupstimLocksublogPowernonLPCCh';

dataPath = [dirPath 'group/Spectral_Data/'];
load([dataPath fileName '.mat'])

% channel groupings/ROIs
nChans = numel(data.ROIid);
groups = zeros(nChans,1);
groups(data.ROIid==1) = 1; %IPS
groups(data.ROIid==2) = 2; %SPL
groups(data.ROIid==3) = 3; %AG

%%
AnalysisNum =3;

% Analysis 1
switch AnalysisNum
    
%% analysis set 1: univariate activity    
    case 1
        % Absstract vs Concrete
        info            =[];
        info.Analysis   = 'stimHGP_AbsConcByROI';
        cond1           = data.conds(:,1); % abstract
        cond2           = data.conds(:,2); % concrete
        info.legend 	= {'Abs','Conc'};
        
        
        % selection for statistics
        X 		= subSelectByCell(data.BinERP,cond1);
        Y 		= subSelectByCell(data.BinERP,cond2);
        
        stats = BiVariateStatsWrapper(X,Y,groups);
        
        info.rownames 	= cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        info.groupNames = {'IPS','SPL','AG'};
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        printStats(stats,info)
        
        % selection for plotting
        addpath Plotting/
        
        X 		= subSelectByCell(data.ERP,cond1);
        Y 		= subSelectByCell(data.ERP,cond2);
        
        t = data.trialTime;
        Z =[];
        Z{1} = cell2mat(cellfun(@mean,X,'uniformoutput',0));
        Z{2} = cell2mat(cellfun(@mean,Y,'uniformoutput',0));
        
        info.cols 		= 'oc';
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.fileName 	= strcat(info.Analysis,'AllChans');
        plotWrapper(Z,t,info)
        
        % by ROI
        for iROI = 1:3
            chans = groups==iROI;
            Z =[];
            Z{1} = cell2mat(cellfun(@mean,X(chans),'uniformoutput',0));
            Z{2} = cell2mat(cellfun(@mean,Y(chans),'uniformoutput',0));
            
            info.fileName 	= sprintf('%s%s',info.Analysis,info.groupNames{iROI});
            plotWrapper(Z,t,info)
        end
        
    case 2
        % Subjective Absstract vs Concrete
        info            =[];
        info.Analysis   = 'stimHGP_RespAbsConcByROI';
        cond1           = data.conds(:,3); % abstract
        cond2           = data.conds(:,4); % concrete
        info.legend 	= {'Resp. Abs','Resp. Conc'};
        
        % selection for statistics
        X 		= subSelectByCell(data.BinERP,cond1);
        Y 		= subSelectByCell(data.BinERP,cond2);
        
        stats = BiVariateStatsWrapper(X,Y,groups);
        
        info.rownames 	= cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        info.groupNames = {'IPS','SPL','AG'};
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        printStats(stats,info)
        
        % selection for plotting
        addpath Plotting/
        
        X 		= subSelectByCell(data.ERP,cond1);
        Y 		= subSelectByCell(data.ERP,cond2);
        
        t = data.trialTime;
        Z =[];
        Z{1} = cell2mat(cellfun(@mean,X,'uniformoutput',0));
        Z{2} = cell2mat(cellfun(@mean,Y,'uniformoutput',0));
        
        info.cols 		= 'oc';
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.fileName 	= strcat(info.Analysis,'AllChans');
        plotWrapper(Z,t,info)
        
        % by ROI
        for iROI = 1:3
            chans = groups==iROI;
            Z =[];
            Z{1} = cell2mat(cellfun(@mean,X(chans),'uniformoutput',0));
            Z{2} = cell2mat(cellfun(@mean,Y(chans),'uniformoutput',0));
            
            info.fileName 	= sprintf('%s%s',info.Analysis,info.groupNames{iROI});
            plotWrapper(Z,t,info)
        end
        
    case 3
        % Correct Absstract vs Concrete
        info            =[];
        info.Analysis   = 'stimHGP_CorrectAbsConcByROI';
        cond1 	= cellfun(@and,data.conds(:,1),data.conds(:,3),'uniformoutput',0); % correct abstract
        cond2 	= cellfun(@and,data.conds(:,2),data.conds(:,4),'uniformoutput',0); % correct concrete
        info.legend 	= {'Correct Abs','Correct Conc'};
        
        % selection for statistics
        X 		= subSelectByCell(data.BinERP,cond1);
        Y 		= subSelectByCell(data.BinERP,cond2);
        
        stats = BiVariateStatsWrapper(X,Y,groups);
        
        info.rownames 	= cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        info.groupNames = {'IPS','SPL','AG'};
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        printStats(stats,info)
        
        % selection for plotting
        addpath Plotting/
        
        X 		= subSelectByCell(data.ERP,cond1);
        Y 		= subSelectByCell(data.ERP,cond2);
        
        t = data.trialTime;
        Z =[];
        Z{1} = cell2mat(cellfun(@mean,X,'uniformoutput',0));
        Z{2} = cell2mat(cellfun(@mean,Y,'uniformoutput',0));
        
        info.cols 		= 'oc';
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.fileName 	= strcat(info.Analysis,'AllChans');
        plotWrapper(Z,t,info)
        
        % by ROI
        for iROI = 1:3
            chans = groups==iROI;
            Z =[];
            Z{1} = cell2mat(cellfun(@mean,X(chans),'uniformoutput',0));
            Z{2} = cell2mat(cellfun(@mean,Y(chans),'uniformoutput',0));
            
            info.fileName 	= sprintf('%s%s',info.Analysis,info.groupNames{iROI});
            plotWrapper(Z,t,info)
        end

%% activity set 2: correlation to test RTs activity
    case 4 
        %% overall correlation
        info            =[];
        info.Analysis   = 'stimHGP_testRTcorr';
        info.groupNames = {'IPS','SPL','AG'};
        info.legend 	= info.groupNames;
        
        X = CorrByCell(data.BinERP,cellfun(@log10,data.testRTs,'uniformoutput',0));
        X = cellfun(@atanh,X,'uniformoutput',0); % convet to continous
        X = [X{:}]';
        
        stats = statsWrapperMat(X,groups);
        info.rownames 	= cellstr(strcat('BinCenter',num2str(mean(data.Bins,2),3)));
        
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Logs/';
        printStats(stats,info)
        
        t = mean(data.Bins,2);
        Z =[];
        Z{1} = X(groups==1,:);
        Z{2} = X(groups==2,:);
        Z{3} = X(groups==3,:);
        
        info.cols 		= 'rbg';
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.fileName 	= strcat(info.Analysis,'AllROIs');
        plotWrapper(Z,t,info)
        
        t = mean(data.Bins,2);
        Z =[];
        Z{1} = repmat(stats.groupScores(1,:),[3,1]); % hard coded for now.
        Z{2} = repmat(stats.groupScores(2,:),[3,1]);
        Z{3} = repmat(stats.groupScores(3,:),[3,1]);
        
        info.cols 		= 'rbg';
        info.savePath 	= '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Results/Plots/';
        info.fileName 	= strcat(info.Analysis,'AllROIsScores');
        plotWrapper(Z,t,info)
            
end

