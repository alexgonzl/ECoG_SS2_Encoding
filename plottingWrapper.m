function  plottingWrapper(figNum)
% plots for ECoG_SS2enconding.

addpath ~/Documents/ECoG_SS2_Encoding/lib/
addpath(genpath('~/Documents/ECoG_SS2_Encoding/lib/contourfcm'))
addpath ~/Documents/ECoG_SS2_Encoding/Analysis/
addpath ~/Documents/ECoG_SS2_Encoding/Plotting/
addpath '~/Documents/ECoG_SS2_Encoding/lib/fdr_bh/'


switch figNum
    case 1
        % I. electrode coverage
        savePath = '~/Google Drive/Research/ECoG_SS2e/Plots/a1/';
        renderChanCortexSS2e(savePath)
    case 2
        % II. Behavioral
        savePath = '~/Google Drive/Research/ECoG_SS2e/Plots/a2/';
        if ~exist(savePath,'dir'), mkdir(savePath), end;
        %behavPlots(savePath);
        behavPlots(savePath);
    case 3
        % III. Activity Plots
        % ROI Encoding Activity
        %lockType     = {'preStim','preStim2','stim','RT'};
        lockType     = {'preStim'};%,'preStim2','stim','RT'};
        for lt = lockType
            opts=[];
            opts.lock = lt{1};
            opts.pThr = 0.05;
            opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a3/' opts.lock '/'];
            if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
            SS2e_MBAnalysis(lt{1},1,opts);
            %SS2e_MBAnalysis(lt{1},3,opts);
            %SS2e_MBAnalysis(lt{1},6,opts);
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
        % GLMs PCA trial Analysis (R^2)
        opts=[];
        opts.rThr = 0.3;
        opts.pThr = 0.01;
        opts.tThr = 1;
        opts.plot1 = 1; opts.plot2 = 0; opts.plot3 = 0;
        opts.plot4 = 0;
        opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        lockType     = {'preStim','preStim2','stim','RT'};
        for lt = lockType
            opts.lock = lt{1};
            opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a6/' opts.lock '/'];
            if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end
            GLM_PCAtrialPlots(opts)
        end
        
    case 7
        % GLMs PCA components analyses
        %         opts=[];
        %         opts.tThr = 1.5;
        %         opts.rThr = 0.15;
        %         opts.pThr = 0.01;
        %         opts.plot1 = 1; opts.plot2 = 0; opts.plot3 = 0;
        %         opts.plot4 = 0; opts.plot5 = 0; opts.plot6 = 0;
        %         opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        %         lockType     = {'preStim2','stim','RT'};
        %
        %         for lt = lockType
        %             opts.lock = lt{1};
        %             opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a7/' opts.lock '/'];
        %             if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
        %             GLM_PCACompPlots(opts)
        %         end
    case 8
        % PCA Components Plots
        opts=[];
        opts.rThr = 0.2;
        opts.pThr = 0.01;
        opts.plot1 = 1; opts.plot2 = 1; opts.plot3 = 0;
        opts.plot4 = 0; opts.plot5 = 1; opts.plot6 = 1;
        opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        lockType     = {'preStim'};%,'preStim2','stim','RT'};
        
        %lockType     = {'RT'};
        %fileName = ['PCATrialDecomp-MBAnalysis2_Kmeans' opts.lock 'sublogPowernonLPCch'];
        
        for lt = lockType
            
            opts.lock = lt{1};
            opts.fileName =['PCATrialDecomp-MBAnalysis3' opts.lock 'sublogPowernonLPCch'];
            opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a8/' opts.lock '/'];
            if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
            
            GLM_PCACompPlots2(opts)
        end
    case 9
        % PCA Components Plots
        opts=[];
        opts.rThr = 0.2;
        opts.pThr = 0.01;
        opts.plot1 = 0; opts.plot2 = 0; opts.plot3 = 0;
        opts.plot4 = 0; opts.plot5 = 1; opts.plot6 = 0;
        opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        lockType     = {'preStim'};%,'preStim2','stim','RT'};
        
        for lt = lockType
            opts.lock = lt{1};
            opts.fileName =['PCATrialDecomp-MBAnalysis_BigBins' opts.lock 'sublogPowernonLPCch'];
            opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a9/' opts.lock '/'];
            if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
            
            GLM_PCACompPlots2(opts)
        end
        
    case 10
        % PCA on single bands
        opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        lockType     = {'preStim2','stim','RT'};
        %lockType     = {'RT'};
        
        for lt = lockType
            opts.lock = lt{1};
            opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a10/' opts.lock '/'];
            if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
            GLM_PCA_SB_CompPlots(opts)
        end
        
    case 11
        % Pre/Post Stim Analyses Plots
        opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        lockType     = {'preStim'};
        
        for lt = lockType
            for ps = 0 %1]
                opts.lock = lt{1};
                opts.poststim500=ps;
                opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a11/' opts.lock '/'];
                if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
                PrePostActRT_Plots(opts)
            end
        end
    case 12
        % colorbar!
        figure(); set(gcf,'paperpositionmode','auto','color','white')
        AR = [80 200];
        set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
        set(gcf,'position',[100,100,AR]); % 100pt all around margin
        ax   = axes('position',[0.15 0.1 0.50 0.8]);
        
        CM=brewermap(9,'*RdBu');
        colorbar2(ax,CM)
        savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a12/'];
        if ~exist(savePath,'dir'), mkdir(savePath), end;
        fN = 'rdb9_colorbar';
        print(gcf,'-dpng',[savePath fN])
        
        cla;
        CM=brewermap(7,'*RdBu');
        colorbar2(ax,CM)
        savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a12/'];
        if ~exist(savePath,'dir'), mkdir(savePath), end;
        fN = 'rdb7_colorbar';
        print(gcf,'-dpng',[savePath fN])
    case 13
        % Pre/Post Stim Analyses Plots
        opts.dataPath = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        lockType     = {'preStim'};
        
        for lt = lockType
            for ps = 0 %1]
                opts.lock = lt{1};
                opts.poststim500=ps;
                opts.savePath = ['~/Google Drive/Research/ECoG_SS2e/Plots/a13/' opts.lock '/'];
                if ~exist(opts.savePath,'dir'), mkdir(opts.savePath), end;
                PrePostActRT_RegionPairs_Plots(opts)
            end
        end
    otherwise
        error('figure has not been implemented.')
end
