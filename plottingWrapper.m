function  plottingWrapper(figNum)
% plots for ECoG_SS2enconding.

addpath ~/Documents/ECoG_SS2_Encoding/lib/
addpath(genpath('~/Documents/ECoG_SS2_Encoding/lib/contourfcm'))
addpath ~/Documents/ECoG_SS2_Encoding/Analysis/
addpath ~/Documents/ECoG_SS2_Encoding/Plotting/

switch figNum
    case 1
        % I. electrode coverage
        savePath = '~/Google Drive/Research/ECoG_SS2e/Plots/a1/';
        renderChanCortexSS2e(savePath)
    case 2
        % II. Behavioral
        savePath = '~/Google Drive/Research/ECoG_SS2e/Plots/a2/';
        if ~exist(savePath,'dir'), mkdir(savePath), end;
        behavPlots(savePath);
    case 3
        % III. Activity Plots
        % ROI Encoding Activity
        lockType     = {'preStim2','stim','RT'};
        for lt = lockType
            opts=[];
            opts.lock = lt{1};
            opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a3/' opts.lock '/'];
            if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
            SS2e_MBAnalysis(lt{1},1,opts)
        end
    case 4
        % IV. Activity Kmeans
        opts.Kmeans.Thtr = 0;
        opts.Kmeans.K    = 3;
        lockType     = {'preStim2','stim','RT'};
        for lt = lockType
            opts.lock = lt{1};
            opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a4/' opts.lock '/'];
            if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
            SS2e_MBAnalysis(lt{1},2,opts)
        end
    case 5
        %V PCA trial analysis
        opts=[];
        opts.rThr = 0.1;
        opts.pThr = 0.05;
        opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        lockType     = {'preStim2','stim','RT'};
        for lt = lockType
            opts.lock = lt{1};
            opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a5/' opts.lock '/'];
            if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
            PCAtrialDecomPlots(opts)
        end
        
    case 6
        % GLMs PCA trial Analysis
        opts=[];
        opts.rThr = 0.3;
        opts.pThr = 0.01;
        opts.tThr = 1;
        opts.plot1 = 0; opts.plot2 = 1; opts.plot3 = 0;
        opts.plot4 = 1;
        opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        lockType     = {'preStim2','stim','RT'};        
        for lt = lockType
            opts.lock = lt{1};
            opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a6/' opts.lock '/'];
            if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end
            GLM_PCAtrialPlots(opts)
        end
        
    case 7
        % GLMs PCA components analyses
        opts=[];
        opts.tThr = 1.6;
        opts.pThr = 0.01;
        opts.plot1 = 0; opts.plot2 = 0; opts.plot3 = 0;
        opts.plot4 = 0; opts.plot5 = 0; opts.plot6 = 1;
        opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        lockType     = {'preStim2','stim','RT'};
        
        for lt = lockType
            opts.lock = lt{1};
            opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a7/' opts.lock '/'];
            if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
            GLM_PCACompPlots(opts)
        end
        
    otherwise
        error('figure has not been implemented.')
end
