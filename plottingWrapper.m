function  plottingWrapper(figNum)
% plots for ECoG_SS2enconding. 

addpath ~/Documents/ECoG_SS2_Encoding/lib/
addpath(genpath('~/Documents/ECoG_SS2_Encoding/lib/contourfcm'))
addpath ~/Documents/ECoG_SS2_Encoding/Analysis/

switch figNum
    case 1
    % I. electrode coverage
    savePath = '~/Google Drive/Research/ECoG_SS2e/Plots/';
    renderChanCortexSS2e(savePath)
	case 2
	% II. Behavioral
	behavPlots();
	case 3
	% III. Activity Plots
	% ROI Encoding Activity 
    lockType     = {'preStim2','stim','RT'};
    for lt = lockType
        SS2e_MBAnalysis(lt{1},1)
    end
    case 4
    % IV. Activity Kmeans
    opts.Kmeans.Thtr = 0;
    opts.Kmeans.K    = 3;    
     lockType     = {'preStim2','stim','RT'};
    for lt = lockType
        SS2e_MBAnalysis(lt{1},2,opts)
    end
    case 5
    %V PCA trial analysis
    opts=[];
    opts.rThr = 0.1;
    opts.pThr = 0.05;
    opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
    opts.savePath = '~/Google Drive/Research/ECoG_SS2e/Plots/PCAtrials/';
    lockType     = {'preStim2','stim','RT'};
    
    for lt = lockType        
        opts.lock = lt{1};
        PCAtrialDecomPlots(opts)
    end
    
    case 6
    % GLMs PCA trial Analysis
    opts=[];
    opts.rThr = 0.3;
    opts.pThr = 0.05;
    opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
    opts.savePath = '~/Google Drive/Research/ECoG_SS2e/Plots/PCAtrials/';
    lockType     = {'preStim2','stim','RT'};
    
    for lt = lockType        
        opts.lock = lt{1};
        GLM_PCAtrialPlots(opts)
    end
    
    

	otherwise
		error('figure has not been implemented.')
end
