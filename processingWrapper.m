%% ecog processing wrapper
% copied from retrieval analyses to encoding analysis
% May 7, 2014

%% preprocess data
addpath PreProcessing/
addpath lib/
addpath behavior/

subjects = {'16b','17b','18','19','24','28','29','30'};
for s = subjects;
    preProcessRawData(s{1})
end

%% re-reference

% updated 2/25/16 to include SMG and TPJ as LPC regions in subjChanInfo  
addpath PreProcessing/
addpath lib/
dataPath ='~/Google Drive/Research/ECoG_SS2e/data/';


subjects = {'16b','17b','18','19','24','28','29','30'};
reference = 'nonLPCCh'; nRefChans = 10; % using this nonLPC for now % 5/7/2014
%reference = 'origCAR'; nRefChans = 0;
%reference = 'nonLPCleasL1TvalCh'; nRefChans = 10;
%reference = 'nLPClowEvokedVar'; nRefChans = 25;

for s = subjects
    fprintf( '\nRe-Referencing Channels for Subject %s \n',s{1});
    dataIn = load([dataPath '/' s{1} '/preProcessed/data.mat']);
    if strcmp(reference,'nLPClowEvokedVar')
        dataIn.data.refInfoFile = [dataPath s{1} '/ERP_Data/ERPsstimLocksubAmporigCAR.mat'];
        data2=load(dataIn.data.refInfoFile);
        dataIn.data.evokedVar =data2.data.evokedVar; clear data2;
    end
    data = reReferenceData(dataIn.data,reference,nRefChans);
    save([dataPath '/' s{1} '/preProcessed/data' reference '.mat'],'data')
    fprintf( 'Re-Referencing for Subject %s Completed. \n',s{1});
end

%% calc ERPs
% addpath PreProcessing/
% addpath Analysis/
% addpath lib/
% dataPath = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
% 
% subjects = {'16b','17b','18','19','24','28','29','30'};
% %reference = 'origCAR';
% %reference = 'nLPClowEvokedVar';
% reference = 'nonLPCCh';
% 
% lockType = {'preStim','stim'};
% 
% 
% analysisType = 'Amp';%{'Amp','Power', 'logPower'};
% baselineType = 'sub';%{'rel','sub'}
% 
% for s = subjects
%     dataIn = load([dataPath s{1} '/preProcessed/data' reference  '.mat']);    
%     dataIn.data.analysisType    = analysisType;
%     dataIn.data.baselineType    = baselineType;
%     for lt = lockType
%         dataIn.data.lockType        = lt{1};
%         switch lt{1}
%             case {'stim','preStim'}
%                 data = calcERP(dataIn.data);
%                 if ~exist([dataPath s{1}  '/ERP_Data/'],'dir'); mkdir([dataPath s{1}  '/ERP_Data/']);end;
%                 save([dataPath s{1}  '/ERP_Data/ERPs' lt{1} baselineType analysisType reference '.mat'],'data')
%             case 'RT' % to do...
%                 % load stim locked data
%                 %             dataIn2 = load([dataPath 'ERP_Data/subj' s{1} '/ERPsstimLock' baselineType analysisType ...
%                 %                 reference num2str(nRefChans) '.mat']);
%                 %             dataIn.data.baseLineMeans = dataIn2.data.baseLineMeans; dataIn2=[];
%                 %             data = calcERP(dataIn.data);
%                 %             save([dataPath 'ERP_Data/subj' s{1} '/ERPsRTLock' baselineType analysisType ...
%                 %                 reference num2str(nRefChans) '.mat'],'data')
%         end
%     end
%     fprintf('calc ERP done for subj %s\n', s{1})
% end

%% bandpass data

addpath PreProcessing/
addpath Analysis/
addpath lib/
dataPath = '~/Google Drive/Research/ECoG_SS2e/data/';

subjects = {'16b','17b','18','19','24','28','29','30'};
%reference = 'nLPClowEvokedVar';
reference = 'nonLPCCh';

for s = subjects
    dataIn = load([dataPath s{1} '/preProcessed/data' reference  '.mat']);
    for band = {'delta','theta','alpha','beta','lgam','hgam'};
        data = dataDecompose(dataIn.data,band{1});
        if ~exist([dataPath s{1}  '/Spectral_Data/continous/'],'dir');
            mkdir([dataPath s{1}  '/Spectral_Data/continous/']);
        end;
        save([dataPath s{1} '/Spectral_Data/continous/BandPass' band{1} reference '.mat'],'data')
    end
end

%% calc ERSPs

addpath PreProcessing/
addpath Analysis/
addpath lib/
dataPath = '~/Google Drive/Research/ECoG_SS2e/data/';

subjects = {'16b','17b','18','19','24','28','29','30'};
reference = 'nonLPCch';
%reference = 'nLPClowEvokedVar';

lockType     = {'preStim2','stim','RT'};

analysisType = 'logPower';%{'Amp','Power', 'logPower'};
baselineType = 'sub';%{'rel','sub'}

for s = subjects
    for band = {'delta','theta','alpha','beta','lgam','hgam'};
        dataIn = load([dataPath s{1} '/Spectral_Data/continous/BandPass' band{1} reference '.mat']);
        dataIn.data.analysisType    = analysisType;
        dataIn.data.baselineType    = baselineType;
            
        for lt = lockType            
            dataIn.data.lockType        = lt{1};
            
            if strcmp(lt{1},'RT')
                % load the stim first
                stimdata = load( [dataPath s{1}  '/Spectral_Data/' band{1} '/ERSPs' band{1} 'stim' baselineType analysisType ...
                    reference '.mat'],'data');
                dataIn.data.baseLineMeans = stimdata.data.baseLineMeans;
            end
            data = calcERSP(dataIn.data);

            if ~exist([dataPath s{1}  '/Spectral_Data/' band{1} '/'],'dir');
                mkdir([dataPath s{1}  '/Spectral_Data/' band{1} '/']);
            end;
            save([dataPath s{1}  '/Spectral_Data/' band{1} '/ERSPs' band{1} lt{1} baselineType analysisType ...
                reference '.mat'],'data')
        end
    end
end

%% group data for individual bands
addpath Analysis/
addpath lib/

bands        = {'delta','theta','alpha','beta','lgam','hgam'};
lockType     = {'preStim2','stim','RT'};
%lockType     = {'RT'};

opts                = [];
opts.hems           = 'all';
%opts.reference      = 'nLPClowEvokedVar';
opts.reference      = 'nonLPCch';
%subjects = {'16b','17b','18','19','24','28','29','30'};
opts.subjects       = {'16b','17b','18','19','24','28','29','30'};
opts.dataPath       = '~/Google Drive/Research/ECoG_SS2e/data/';

for lt = lockType
    opts.lockType       = lt{1};
    for ba = 1:numel(bands)
        opts.band   = bands{ba};
        data        = groupLPCData(opts);
        fileName    = [opts.hems data.prefix 'Group' data.extension];
        
        if strcmp(opts.band,'erp')            
            savePath = [opts.dataPath 'group/ERP_Data/'];
            savePath2 = '~/Google Drive/Research/ECoG_SS2e/data_results/';
        else
            savePath = [opts.dataPath 'group/Spectral_Data/'];
            savePath2 = ['~/Google Drive/Research/ECoG_SS2e/data_results/' opts.lockType '/'];
        end
        
        if ~exist(savePath2,'dir'); mkdir(savePath2); end;
        
        save([savePath fileName '.mat'],'data')
        save([savePath2 fileName '.mat'],'data')
        fprintf('grouping data completed for %s\n',opts.band)
    end
end


%% group across frequency bands

addpath Analysis/
addpath lib/

bands        = {'delta','theta','alpha','beta','lgam','hgam'};
lockType     = {'preStim2','stim','RT'};
%lockType     = {'preStim2'};
opts                = [];
opts.hems           = 'all';
opts.bands          = bands;
%opts.reference      = 'nLPClowEvokedVar';
opts.reference      = 'nonLPCch';

for lt = lockType
    opts.lockType       = lt{1};
    opts.dataPath       = ['~/Google Drive/Research/ECoG_SS2e/data_results/' opts.lockType '/'];
    data        = groupLPCDataMultiBand(opts);        
    fileName   = [opts.hems 'MBAnalysis' data.extension];
    save([opts.dataPath fileName '.mat'],'data')
    fprintf('grouping data completed for %s\n',lt{1})    
end

%% obtain channel locations & plot locs

dataSS2ret=load('~/Google Drive/Research/ECoG Manuscript/data/allERSPshgamGroupstimLocksublogPowernonLPCleasL1TvalCh10.mat');
dataSS2enc=load('~/Google Drive/Research/ECoG_SS2e/data_results/stim/allERSPsdeltaGroupstimsublogPowernonLPCch.mat');
SS2e_subjects       = {'16b','17b','18','19','24','28','29','30'};

elecLocs            = [];
elecLocs.MNILocs    = [];
elecLocs.MNIcortex  = dataSS2ret.data.MNIcortex;
elecLocs.lMNIcortex  = dataSS2ret.data.lMNIcortex;
elecLocs.rMNIcortex  = dataSS2ret.data.rMNIcortex;
elecLocs.ROIid      = dataSS2enc.data.ROIid;
elecLocs.SubjChans   =dataSS2enc.data.subjChans;
elecLocs.hemChans    =dataSS2enc.data.hemChan;
for ss = 1:numel(SS2e_subjects)
    %SS2retSubjID = find(strcmp(dataSS2ret.data.options.subjects,SS2e_subjects{ss}));
    load(['lib/elecLocs/subj' SS2e_subjects{ss} '_mni_elcoord_corrected.mat'],'mni_elcoord');
    chans=mni_elcoord(dataSS2enc.data.LPCchanId(elecLocs.SubjChans==ss),:);    
    elecLocs.MNILocs     = [elecLocs.MNILocs;chans];     
end
dataPath       = '~/Google Drive/Research/ECoG_SS2e/data_results/Renderings/';
fileName       = 'electrodeLocs';
save([dataPath fileName '.mat'],'elecLocs');
%% Kmeans temporo-spectral analyses

% addpath Analysis
% addpath lib
% 
% lockType     = {'preStim2','stim','RT'};
% %analysis     = {'activity','studyRT','testRT'};
% analysis     = {'activity'};
% 
% opts = [];
% opts.TvalChanThr    = 0; % # inclusion criteria for channels; zero includes all chans
% opts.numClusters    = 3; % K opa
% opts.replicates     = 100;
% opts.reference      = 'nonLPCch';
% opts.dataPath = ['~/Google Drive/Research/ECoG_SS2e/data_results/'];
% 
% for lt = lockType
%     for an = analysis
%         opts.lockType = lt{1};
%         
%         opts.analysis = an{1};        
%         extension = [opts.lockType 'sublogPower' opts.reference];
%         fileName =  ['allMBAnalysis' extension];
%     
%         load([opts.dataPath opts.lockType '/' fileName '.mat'])
%         opts.chans  = sum(data.chBinScore>opts.TvalChanThr)>1;
%         out = kmeansAnalyses(data,opts);
%     
%         save([opts.dataPath 'Kmeans/' opts.lockType '_' opts.analysis '_Tthr' ...
%             strrep(num2str(opts.TvalChanThr),'.','p') '_Kmeans' num2str(opts.numClusters)...
%             'MultiBand' extension '.mat'],'out')     
%     end
% end 
%% PCA trial analysis
addpath Analysis/
addpath lib/

%lockType     = {'preStim2','stim','RT'};
lockType     = {'preStim2','stim','RT'};

opts                = [];
opts.hems           = 'all';
opts.nComps         = 12;
opts.reference      = 'nonLPCch';
opts.dataPath       = '~/Google Drive/Research/ECoG_SS2e/data_results/';
for lt = lockType
    opts.lockType       = lt{1};
    extension           = [opts.lockType 'sublogPower' opts.reference];
    fileName            = [opts.hems 'MBAnalysis' extension];
    load([opts.dataPath opts.lockType '/' fileName '.mat'])
    
    out     = PCATrialDecomp(data,opts);
    
    fileName            = ['PCATrialDecomp-MBAnalysis' extension]; 
    save([opts.dataPath opts.lockType '/' fileName '.mat'],'out')
    fprintf('PCA trial decomp completed for %s\n',lt{1})    
end

%% multiband ITC

addpath PreProcessing/
addpath Analysis/
addpath lib/
addpath lib/CircStats/
dataPath = '~/Google Drive/Research/ECoG_SS2e/data/';

subjects = {'16b','17b','18','19','24','28','29','30'};
reference = 'nonLPCch';
lockType     = {'preStim2','stim','RT'};

opts = []; 
opts.nFreq = 30;
opts.downsample = 2;
opts.NumComponents = 12;
for s = subjects    
    dataIn = load([dataPath s{1} '/preProcessed/data' reference  '.mat']);            
    for lt = lockType            
        opts.lockType        = lt{1};            
        data = multiBandITC(dataIn.data,opts);

        if ~exist([dataPath s{1}  '/ITC_Data/'],'dir');
            mkdir([dataPath s{1}  '/ITC_Data/']);
        end;
        save([dataPath s{1}  '/ITC_Data/multiBandITC_N' num2str(opts.nFreq) '_' ...
            lt{1} reference '.mat'],'data')                
        disp(sprintf('processing completed for %s, subj %s',lt{1},s{1}));
        data = multiBandITC_GLMPCA(data,opts);
        save([dataPath s{1}  '/ITC_Data/multiBandITC_GLMPCA_N' num2str(opts.nFreq) '_' ...
            lt{1} reference '.mat'],'data')   
        disp(sprintf('Phase PCA-GLM completed for %s, subj %s',lt{1},s{1}));
    end      
end

% modulation index computation using multi bands
addpath PreProcessing/
addpath Analysis/
addpath lib/
addpath lib/CircStats/
dataPath = '~/Google Drive/Research/ECoG_SS2e/data/';

subjects = {'16b','17b','18','19','24','28','29','30'};
reference = 'nonLPCch';
lockType     = {'preStim2','stim','RT'};

for s = subjects    
    
    for lt = lockType            
        opts.lockType        = lt{1};            
        ampdata = load([dataPath s{1} '/Spectral_Data/hgam/ERSPshgam' lt{1} 'sublogPower' reference '.mat' ]);                    
        ampdata = ampdata.data;
        phdata  =load([dataPath s{1} '/ITC_Data/multiBandITC_N30_' lt{1} reference '.mat']);                    
        phdata = phdata.data;
        if ~exist([dataPath s{1}  '/MI_Data/'],'dir');
            mkdir([dataPath s{1}  '/MI_Data/']);
        end;
        data = calcMI(ampdata,phdata);
        save([dataPath s{1}  '/MI_Data/modIndex_' lt{1} reference '.mat'],'data')                
        disp(sprintf('processing completed for %s, subj %s',lt{1},s{1}));
        
    end      
end

%% multiband ITC PCA-GLM



% group subjects.

% opts = []; 
% opts.nFreq = 30;
% 
% for lt = lockType            
%     opts.lockType = lt{1};
%     data = groupMultiBandITCData(opts);
% end
%% lasso and ridge analysis

%%
%  Plot Correlation time-courses per band and ROI
% 
% addpath Analysis/
% addpath lib/
% addpath Plotting/
% 
% bands        = {'hgam','theta','alpha'};
% lockTypes     = {'preStim','stim','RT'};
% hemispheres  = {'left','right','all'};
% rois            = {'IPS','SPL'};
% cols            = {'r','b'};
% dataPath       = '~/Google Drive/Research/ECoG_SS2e/data_results/';
% 
% info           = [];
% info.alpha      = 0.001;
% info.savePath  = '~/Google Drive/Research/ECoG_SS2e/Plots/';
% 
% for lt = 1:numel(lockTypes)
%     lockType       = lockTypes{lt};
%     if strcmp(lockType,'RT')
%         info.yAxisRightLoc=1;
%     else
%         info.yAxisRightLoc=0;
%     end
%     for ba = 1:numel(bands)
%         band   = bands{ba};
%         fileName    = ['allERSPs' band 'Group' lockType 'sublogPowernLPClowEvokedVar'];
%         load([dataPath fileName '.mat'],'data')
%         b = mean(data.Bins,2);
%         t = data.trialTime;
%         for hm = 1:numel(hemispheres)
%             hem = hemispheres{hm};
%             switch hem
%                 case 'left'
%                     ROIch{1} = data.hemChan==1 & data.ROIid==1;
%                     ROIch{2} = data.hemChan==1 & data.ROIid==2;
%                 case 'right'
%                     ROIch{1} = data.hemChan==2 & data.ROIid==1;
%                     ROIch{2} = data.hemChan==2 & data.ROIid==2;
%                     
%                 case 'all'
%                     ROIch{1} = data.ROIid==1;
%                     ROIch{2} = data.ROIid==2;
%             end
%             for rr = 1:numel(rois)
%                 info.cols = cols{rr};
%                 % mean activity
%                 X = [];
%                 X{1} = data.meanChResp(ROIch{rr},:);
%                 info.fileName = [hem '_' band '_' lockType '_' rois{rr} '_meanChResp'];
%                 plotWrapper(X,t,info)
%                 
%                 % correlation of activity to semantic desicion
%                 X{1} = data.dataToStudyRTsCorr(ROIch{rr},:);
%                 [~,info.PVals] = ttest(X{1});
%                 info.fileName = [hem '_' band '_' lockType '_' rois{rr} '_CorrToSemanticRT'];
%                 plotWrapper(X,b,info)
%                 
%                 % correlation of activity to retrieval RT
%                 X{1} = data.dataToTestRTsCorr(ROIch{rr},:);
%                 [~,info.PVals] = ttest(X{1});
%                 info.fileName = [hem '_' band '_' lockType '_' rois{rr} '_CorrToRetrievalRT'];
%                 plotWrapper(X,b,info)
%             end
%         end
%     end
% end

%renderChanCortexSS2e(dataPath);

%
% calc ITC

% addpath PreProcessing/
% addpath Analysis/
% addpath lib/
% dataPath = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
% 
% subjects = {'16b','17b','18','24','28','29','30'};
% %reference = 'nonLPCch';
% reference = 'nLPClowEvokedVar';
% 
% lockType     = {'preStim','stim'};
% 
% analysisType = 'logPower';%{'Amp','Power', 'logPower'};
% baselineType = 'sub';%{'rel','sub'}
% 
% for s = subjects
%     for band = {'delta','theta','alpha','beta','lgam','hgam'};
%         for lt = lockType
%             dataIn = load([dataPath s{1} '/Spectral_Data/continous/BandPass' band{1} reference '.mat']);
%             dataIn.data.lockType        = lt{1};
%             dataIn.data.analysisType    = analysisType;
%             dataIn.data.baselineType    = baselineType;
%             
%            
%             data = calcITC(dataIn.data);
%             if ~exist([dataPath s{1}  '/ITC_Data/' band{1} '/'],'dir');
%                 mkdir([dataPath s{1}  '/ITC_Data/' band{1} '/']);
%             end;
%             save([dataPath s{1}  '/ITC_Data/' band{1} '/ITC' band{1} lt{1} baselineType analysisType ...
%                 reference '.mat'],'data')
%         end
%     end
% end
% 
% 
% %% group ITC
% 
% addpath Analysis/
% addpath lib/
% 
% bands        = {'delta','theta','alpha'};
% lockType     = {'preStim','stim'};%{'preStim','stim','RT'};
% 
% opts                = [];
% opts.hems           = 'all';
% opts.reference      = 'nLPClowEvokedVar';
% %opts.reference      = 'nonLPCch';
% opts.subjects       = {'16b','17b','18','24','28','29','30'};
% opts.dataPath       = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
% 
% for lt = lockType
%     opts.lockType       = lt{1};
%     for ba = 1:numel(bands)
%         opts.band   = bands{ba};
%         data        = groupLPC_ITCData(opts);
%         fileName    = [opts.hems data.prefix 'Group' data.extension];
%         
%         if strcmp(opts.band,'erp')            
%             savePath1 = [opts.dataPath 'group/ERP_Data/'];
%             savePath2 = '~/Google Drive/Research/ECoG_SS2e/data_results/';
%         else
%             savePath1 = [opts.dataPath 'group/ITC_Data/'];
%             savePath2 = '~/Google Drive/Research/ECoG_SS2e/data_results/';
%         end
%         
%         if ~exist(savePath1,'dir'); mkdir(savePath1); end;
%         
%         save([savePath1 fileName '.mat'],'data')
%         save([savePath2 fileName '.mat'],'data')
%         fprintf('grouping data completed for %s\n',opts.band)
%     end
% end
% 
% %% calc MI
% 
% %dateStr = '27-May-2013';
% %subjects = {'16b','18','24','28'};
% dateStr = '17-Jun-2013';
% subjects = {'17b','19','29'};
% reference = 'nonLPCleasL1TvalCh'; nRefChans = 10;
% lockType = 'stim';
% dataPath = '../Results/';
% 
% for s = subjects
%     for AmpBand = {'hgam'}
%         data1 = load([dataPath 'Spectral_Data/subj' s{1} '/BandPassedSignals/BandPass' ...
%             AmpBand{1}  reference num2str(nRefChans) dateStr '.mat']);
%         
%         for PhaseBand = {'delta','theta','alpha','beta'}
%             data2 = load([dataPath 'Spectral_Data/subj' s{1} '/BandPassedSignals/BandPass' ...
%                 PhaseBand{1}  reference num2str(nRefChans) dateStr '.mat']);
%             
%             data = calcMI(data1.data,data2.data,lockType);
%             savePath = [dataPath 'MI_Data/subj' s{1} '/'];
%             
%             if ~exist(savePath,'dir'), mkdir(savePath), end;
%             save([savePath '/MI_AMP' AmpBand{1} '_PHASE' PhaseBand{1} ...
%                 lockType 'Lock' reference num2str(nRefChans) '.mat'],'data')
%         end
%     end
% end
% 
% %% group MI
% 
% % to do ... (noted on Oct. 22 2013 )
% 
% opts = [];
% opts.hems       = 'all';
% opts.lockType = 'stim';
% opts.reference = 'nonLPCleasL1TvalCh'; opts.nRefChans = 10;
% opts.type = 'MI';
% opts.ampBand = 'hgam';
% 
% opts.subjects   = {'16b','18','24','28','17b','19', '29'};
% opts.hemId      = {'l'  ,'l' ,'l' ,'l' ,'r'  ,'r' , 'r'};
% 
% dataPath = '../Results/MI_Data/group/';
% for bands = {'theta'}
%     
%     opts.phaseBand = bands{1};
%     
%     data = groupLPC_MIData(opts);
%     
%     fileName = [opts.hems 'Group' data.fileParams];
%     save([dataPath fileName '.mat'],'data')
%     
% end
% 
% %% plot roi level data
% 
% opts                = [];
% opts.hems           = 'l';
% opts.lockType       = 'RT';
% opts.reference      = 'nonLPCleasL1TvalCh';
% opts.nRefChans      = 10;
% opts.type           = 'power';
% opts.band           = 'hgam';
% opts.smoother       = 'loess';
% opts.acrossWhat     = 'Subjects';
% opts.smootherSpan   = 0.15;
% opts.yLimits        = [-0.6 1.5];
% opts.aRatio         = [500 300];
% 
% opts.subjects       = {'16b','18','24','28','17b','19', '29'};
% opts.hemId          = {'l'  ,'l' ,'l' ,'l' ,'r'  ,'r' , 'r'};
% 
% opts.mainPath = '../Results/' ;
% if strcmp(opts.type,'erp')
%     opts.measType       = 'm';      % {'m','z','c','Zc'}
%     opts.comparisonType = 'ZStat';  % {ZStat,ZcStat}
%     opts.baselineType   = 'sub';
%     opts.analysisType   = 'Amp';
%     opts.dataPath       = [opts.mainPath 'ERP_Data/group/'];
%     opts.preFix         = 'ERPs' ;
%     opts.plotPath       = [opts.mainPath 'Plots/ERPs/'  opts.hems '/'];
%     opts.band           = '';
% elseif strcmp(opts.type,'power')
%     opts.measType       = 'm';     % {'m','z','c','Zc'}
%     opts.comparisonType = 'ZStat'; % {ZStat,ZcStat}
%     opts.baselineType   = 'sub';
%     opts.analysisType   = 'logPower';
%     opts.dataPath       = [opts.mainPath 'Spectral_Data/group/'];
%     opts.preFix         = ['ERSPs' opts.band];
%     opts.plotPath       = [opts.mainPath 'Plots/Spectral/' opts.hems '/'];
% elseif strcmp(opts.type,'ITC')
%     opts.measType       = 'ITC_';
%     opts.comparisonType = 'ITC_Z';
%     opts.baselineType   = '';
%     opts.analysisType   = '';
%     opts.dataPath       = [opts.mainPath 'ITC_Data/group/'];
%     opts.preFix         = ['ITC' opts.band];
%     opts.plotPath       = [opts.mainPath 'Plots/ITC/' opts.hems '/'];
% end
% 
% opts.extension  = [opts.lockType 'Lock' opts.baselineType opts.analysisType opts.reference ...
%     num2str(opts.nRefChans)] ;
% fileName        = ['all' opts.preFix 'Group' opts.extension '.mat'];
% load([opts.dataPath fileName])
% close all
% 
% plotROI_ERPs(data,opts)
% 
% switch opts.lockType
%     case 'RT'
%         opts.timeLims   = [-0.6 0.1];
%         opts.timeStr     = 'n600msTo100ms';
%     case 'stim'
%         opts.timeLims   = [0 1];
%         opts.timeStr     = '0msTo1000ms';
%end
% 
% % opts.aspectRatio = [50 300];
% >>>>>>> 2e86b6b2bcae5544f60700fe62990fe64e1255cf
% % opts.hem        = 1;
% 
% opts.ROInums    = [1];
% plotROI_TC_multiband(opts)
% 
% opts.ROInums    = [2];
% plotROI_TC_multiband(opts)
% 
% opts.yLimits    = [-2.5 1];
% opts.yTicks     = [-2 -1 0];
% conditionBarMultiBandWrapper(opts)
% 
% 
% opts.ROInums    = [2];
% % opts.yLimits    = [-5.5 1];
% % opts.yTicks     = [-4 -2 0];
% opts.yLimits    = [-2.5 1];
% opts.yTicks     = [-2 -1 0];
% conditionBarMultiBandWrapper(opts)
% 
% %% mni group renderings
% 
% addpath Plotting/
% close all
% opts                = [];
% opts.hem            = 'l';
% opts.lockType       = 'RT';
% opts.type           = 'power';
% opts.band           = 'hgam';
% opts.comparisonType = 'BinZStat'; %{'BinTStat','BinpValT','BinpVal','BinZStat','BinpValz','BincondDiff'};
% opts.timeType       = 'Bins';
% opts.reference      = 'nonLPCleasL1TvalCh';
% opts.nRefChans      = 10;
% opts.renderType     = 'SmoothCh';%{'SmoothCh','UnSmoothCh', 'SigChans','SignChans'};
% opts.limitDw        = -4;
% opts.limitUp        = 4;
% opts.absLevel       = 1;
% opts.resolution     = 400;
% opts.avgBins        = [];
% 
% opts.subjects       = {'16b','18','24','28','17b','19', '29'};
% opts.hemId          = {'l'  ,'l' ,'l' ,'l' ,'r'  ,'r' , 'r'};
% 
% opts.mainPath = '../Results/' ;
% if strcmp(opts.type,'erp')
%     opts.baselineType   = 'sub';
%     opts.analysisType   = 'Amp';
%     opts.dataPath       = [opts.mainPath 'ERP_Data/group/'];
%     opts.preFix         = 'ERPs' ;
%     opts.renderPath = [opts.mainPath 'Plots/Renderings/ERPs/'];
% elseif strcmp(opts.type,'power')
%     opts.baselineType   = 'sub';
%     opts.analysisType   = 'logPower';
%     opts.dataPath       = [opts.mainPath 'Spectral_Data/group/'];
%     opts.preFix         = ['ERSPs' opts.band];
%     opts.renderPath     = [opts.mainPath 'Plots/Renderings/Spectral/' opts.band '/'];
% elseif strcmp(opts.type,'ITC')
%     opts.baselineType   = '';
%     opts.analysisType   = '';
%     opts.dataPath       = [opts.mainPath 'ITC_Data/group/'];
%     opts.preFix         = ['ITC' opts.band];
%     opts.renderPath     = [opts.mainPath 'Plots/Renderings/ITC/' opts.band '/'];
% end
% 
% opts.extension = [opts.lockType 'Lock' opts.baselineType opts.analysisType opts.reference ...
%     num2str(opts.nRefChans)] ;
% fileName    = ['all' opts.preFix 'Group' opts.extension '.mat'];
% 
% load([opts.dataPath fileName])
% if ~exist(opts.renderPath,'dir'),mkdir(opts.renderPath),end
% renderERPs(data,opts)
% 
% %% cluster channels
% addpath Plotting/
% addpath Analysis//
% close all
% 
% opts                = [];
% opts.hems            = 'l'; opts.hemNum=1;
% opts.ROIs           = [1 2];
% opts.nClusters      = 4;
% opts.lockType       = 'stim';
% opts.type           = 'power';
% opts.band           = 'hgam';
% opts.dtype          = 'ZStat';
% opts.smoothData     = true;
% opts.reference      = 'nonLPCleasL1TvalCh';
% opts.nRefChans      = 10;
% opts.plotting       = false;
% opts.findSubCluster = true; % only wors for K=2
% opts.resolution     = 400;
% opts.aRatio         = [500 300];
% 
% 
% opts.mainPath = '../Results/' ;
% if strcmp(opts.type,'erp')
%     opts.measType       = 'm';      % {'m','z','c','Zc'}
%     opts.comparisonType = 'ZStat';  % {ZStat,ZcStat}
%     opts.baselineType   = 'sub';
%     opts.analysisType   = 'Amp';
%     opts.dataPath       = [opts.mainPath 'ERP_Data/group/'];
%     opts.preFix         = 'ERPs' ;
%     opts.plotPath       = [opts.mainPath 'Plots/ERPs/'  opts.hems '/'];
%     opts.band           = '';
% elseif strcmp(opts.type,'power')
%     opts.measType       = 'm';     % {'m','z','c','Zc'}
%     opts.comparisonType = 'ZStat'; % {ZStat,ZcStat}
%     opts.baselineType   = 'sub';
%     opts.analysisType   = 'logPower';
%     opts.dataPath       = [opts.mainPath 'Spectral_Data/group/'];
%     opts.preFix         = ['ERSPs' opts.band];
%     opts.plotPath       = [opts.mainPath 'Plots/Spectral/' opts.hems '/'];
% elseif strcmp(opts.type,'ITC')
%     opts.measType       = 'ITC_';
%     opts.comparisonType = 'ITC_Z';
%     opts.baselineType   = '';
%     opts.analysisType   = '';
%     opts.dataPath       = [opts.mainPath 'ITC_Data/group/'];
%     opts.preFix         = ['ITC' opts.band];
%     opts.plotPath       = [opts.mainPath 'Plots/ITC/' opts.hems '/'];
% end
% 
% opts.extension = [opts.lockType 'Lock' opts.baselineType opts.analysisType opts.reference ...
%     num2str(opts.nRefChans)] ;
% fileName    = ['all' opts.preFix 'Group' opts.extension '.mat'];
% 
% load([opts.dataPath fileName])
% out=clusterWrapper(data, opts);
% 
% fileName2 = ['K' num2str(opts.nClusters) 'Clusters' fileName];
% opts.savePath = [opts.dataPath 'clusters/'];
% if ~exist(opts.savePath,'dir'), mkdir(opts.savePath),end
% 
% save([opts.savePath fileName2],'out')
% 
% %% output group data to a csv file and compute stats
% 
% addpath Analysis/
% clc
% 
% opts                = [];
% opts.hems           = 'all';
% opts.lockType       = 'RT';
% opts.reference      = 'nonLPCleasL1TvalCh'; opts.nRefChans = 10;
% opts.type           = 'power'; opts.band = 'hgam';
% %opts.type           = 'power'; opts.band = 'hgam';
% opts.bin            = 'Bin'; % options are{'BigBin', 'Bin'};
% opts.byBlockFlag    = 0;
% 
% switch opts.lockType
%     case 'RT'
%         opts.time   = [-1 0.2];
%         timeStr     = 'n1000msTo200ms';
%     case 'stim'
%         opts.time   = [0 1];
%         timeStr     = '0msTo1000ms';
% end
% 
% 
% dataPath = '../Results/';
% if strcmp(opts.type,'erp')
%     dataPath = [dataPath 'ERP_Data/group/'];
%     fileName = [opts.hems 'ERPsGroup' opts.lockType 'LocksubAmp' ...
%         opts.reference num2str(opts.nRefChans)];
% elseif strcmp(opts.type,'power')
%     dataPath = [dataPath 'Spectral_Data/group/'];
%     fileName = [opts.hems 'ERSPs' opts.band 'Group' opts.lockType  ...
%         'LocksublogPower' opts.reference num2str(opts.nRefChans)];
% elseif strcmp(opts.type,'itc')
%     dataPath = [dataPath 'ITC_Data/group/'];
%     fileName = [opts.hems 'ITC' opts.band 'Group' opts.lockType  ...
%         'Lock' opts.reference num2str(opts.nRefChans)];
% end
% 
% load([dataPath fileName])
% clc
% 
% printStats(data,opts)
% 
% savePath = '../Results/Rdata/';
% blockStr = {'bySubj','byBlock'}; blockStr = blockStr{opts.byBlockFlag+1};
% out = exportLPCData2R(data,opts);
% csvwrite([savePath opts.bin timeStr fileName blockStr '.csv'],single(out));
% 
% %% %%%%%%%%%
% 
% %% decoding
% 
% clear all; close all;
% addpath Classification/
% addpath lib/
% 
% opts                = [];
% opts.lockType       = 'RT';
% opts.reference      = 'nonLPCleasL1TvalCh'; opts.nRefChans = 10;
% %opts.dataType       = 'erp'; opts.bands          = {''};
% opts.dataType       = 'power'; opts.bands          = {'hgam'};%{'delta','theta','alpha','beta','lgam','hgam'};
% %opts.dataType       = 'power'; opts.bands          = {'erp','delta','theta','alpha','beta','lgam','hgam'};
% opts.toolboxNum     = 1;
% 
% % feauture settings
% % timeType distinguishes between decoding individual samples, or taking
% % time  bins
% % options are {'','Bin'};
% opts.timeType       = 'Bin';
% % channelGroupingType makes the distinction between decoding between channels, rois or
% % all channels within LPC
% % options are {'channel','ROI','IPS-SPL','all'};
% opts.channelGroupingType      = 'channel';
% % timeFeatures distinguishes between taking the whole trial or taking a
% % window of time
% % options are {'window','trial'};
% opts.timeFeatures   = 'trial';
% 
% switch opts.lockType
%     case 'RT'
%         opts.timeLims   = [-0.8 0.2];
%         opts.timeStr     = 'n800msTo200ms';
%         %         opts.timeLims   = [-0.6 0.1];
%         %         opts.timeStr     = 'n600msTo100ms';
%     case 'stim'
%         %opts.timeLims   = [-0.2 1];
%         %opts.timeStr     = 'n200msTo1000ms';
%         opts.timeLims   = [0 1];
%         opts.timeStr     = '0msTo1000ms';
% end
% 
% S = ClassificationWrapper(opts);
% 
% savePath = ['../Results/Classification/group/' opts.dataType ...
%     '/' opts.channelGroupingType '/'];
% if ~exist(savePath,'dir'),mkdir(savePath), end
% 
% fileName = ['allSubjsClassXVB' opts.lockType 'Lock' opts.timeStr opts.dataType cell2mat(opts.bands) '_tF' opts.timeFeatures '_tT' ...
%     opts.timeType '_gT' opts.channelGroupingType '_Solver' S.extStr];
% save([savePath fileName],'S')
% 
% %% summarize performance
% 
% addpath Classification/
% 
% opts                = [];
% opts.lockType       = 'RT';
% opts.reference      = 'nonLPCleasL1TvalCh'; opts.nRefChans = 10;
% %opts.dataType       = 'power'; opts.bands          = {'delta','theta','alpha'};
% opts.dataType       = 'power'; opts.bands          = {'hgam'};
% %opts.dataType       = 'power'; opts.bands          = {'erp','hgam'};
% %opts.dataType       = 'power'; opts.bands          = {'delta','theta','alpha','beta','lgam','hgam'};
% opts.toolboxNum     = 1;
% opts.timeType       = 'Bin';
% opts.channelGroupingType      = 'channel';
% opts.timeFeatures   = 'trial';
% opts.extStr         = 'liblinearS0';%'NNDTWK5';%
% 
% switch opts.lockType
%     case 'RT'
%         opts.timeLims   = [-0.8 0.2];
%         opts.timeStr     = 'n800msTo200ms';
%         %         opts.timeLims   = [-0.6 0.1];
%         %         opts.timeStr     = 'n600msTo100ms';
%     case 'stim'
%         %opts.timeLims   = [-0.2 1];
%         %opts.timeStr     = 'n200msTo1000ms';
%         opts.timeLims   = [0 1];
%         opts.timeStr     = '0msTo1000ms';
% end
% 
% dataPath = ['../Results/Classification/group/' opts.dataType ...
%     '/' opts.channelGroupingType '/'];
% 
% fileName = ['allSubjsClassXVB' opts.lockType 'Lock' opts.timeStr opts.dataType cell2mat(opts.bands) '_tF' opts.timeFeatures '_tT' ...
%     opts.timeType '_gT' opts.channelGroupingType '_Solver' opts.extStr];
% 
% load([dataPath fileName])
% S = SummaryClassification(S,opts);
% fileName = ['Sum' fileName];
% % save([dataPath fileName],'S')
% 
% %% plot decoding results
% 
% addpath lib/
% 
% opts                = [];
% opts.lockType       = 'stim';
% opts.scoreType      = 'mBAC'; % RTsLogitCorr mBAC
% opts.accPlots       = true;
% opts.weigthsPlots   = false;
% opts.renderPlot     = true;
% opts.RTcorrPlots    = false;
% opts.stats          = false;
% opts.baseLineY      = 0;
% opts.rendLimits     = [-0.15 0.15];
% opts.resolution     = 400;
% opts.reference      = 'nonLPCleasL1TvalCh'; opts.nRefChans = 10;
% opts.toolboxNum     = 1;
% opts.dataType       = 'power'; opts.bands          = {'hgam'};
% %opts.dataType       = 'power'; opts.bands          = {'delta','theta','alpha'};
% %opts.dataType       = 'power'; opts.bands          = {'theta'};
% %opts.dataType       = 'power'; opts.bands          = {'delta','theta','alpha','beta','lgam','hgam'};
% 
% 
% % accuracy plots options:
% opts.accPlotsOpts.indPoints   = false; % plot individual points
% opts.accPlotsOpts.horizontal  = false; % horizontal plot
% opts.accPlotsOpts.yLimits     = [0.49 0.61]; % horizontal plot
% opts.accPlotsOpts.aspectRatio = [200 600];
% 
% opts.timeType       = 'Bin';
% opts.channelGroupingType      = 'channel';
% opts.timeFeatures   = 'trial';
% opts.extStr         = 'liblinearS0';%'NNDTW_K5';
% 
% dataPath = ['../Results/Classification/group/' opts.dataType ...
%     '/' opts.channelGroupingType '/'];
% 
% fileName = ['SumallSubjsClass' opts.lockType 'Lock' opts.dataType cell2mat(opts.bands) '_tF' opts.timeFeatures '_tT' ...
%     opts.timeType '_gT' opts.channelGroupingType '_Solver' opts.extStr];
% 
% opts.fileName = fileName;
% load([dataPath fileName])
% 
% close all
% plotDecodingAcc(S,opts)
% %plotFA_MISS_relationToACC
% 
% %% plot relationship between channel accuracies
% 
% close all;
% opts                = [];
% opts.lockType1       = 'stim';
% opts.dataType1       = 'power'; opts.bands1        = {'theta','hgam'};
% %opts.dataType1       = 'power'; opts.bands1          = {'delta','theta','alpha','beta','lgam','hgam'};
% opts.lockType2       = 'stim';
% opts.dataType2       = 'power'; opts.bands2        = {'hgam'};
% %opts.dataType2       = 'power'; opts.bands2          = {'delta','theta','alpha','beta','lgam','hgam'};
% 
% opts.subjects       = [1:4]; % left subjects
% opts.ROIs           = [1 2]; % roi 1 and 2, IPS and SPL
% opts.ROIids         = true;  % plot roi colors
% opts.lims = [0.45 0.80];
% 
% opts.timeType       = 'Bin';
% opts.channelGroupingType      = 'channel';
% opts.timeFeatures   = 'trial';
% opts.extStr         = 'liblinearS0';%'NNDTW_K5';
% 
% dataPath = ['../Results/Classification/group/'];
% 
% fileName1 = [ opts.dataType1 '/' opts.channelGroupingType '/' 'SumallSubjsClass' ...
%     opts.lockType1 'Lock' opts.dataType1 cell2mat(opts.bands1) '_tF' opts.timeFeatures ...
%     '_tT' opts.timeType '_gT' opts.channelGroupingType '_Solver' opts.extStr];
% 
% fileName2 = [ opts.dataType2 '/' opts.channelGroupingType '/' 'SumallSubjsClass' ...
%     opts.lockType2 'Lock' opts.dataType2 cell2mat(opts.bands2) '_tF' opts.timeFeatures ...
%     '_tT' opts.timeType '_gT' opts.channelGroupingType '_Solver' opts.extStr];
% 
% opts.fileName = fileName1;
% opts.fileName = fileName2;
% load([dataPath fileName1])
% data1 = S;
% load([dataPath fileName2])
% data2 = S;
% 
% opts.savePath = '/Users/alexg8/Google Drive/Research/ECoG Manuscript/ECoG Manuscript Figures/individualPlotsPDFs';
% close all
% plotACCRelationshipWrapper(data1,data2,opts)
% 
% %% find IPS / SPL time segment matches
% 
% clearvars
% 
% opts                = [];
% opts.lockType       = 'stim';
% opts.reference      = 'nonLPCleasL1TvalCh'; opts.nRefChans = 10;
% %opts.dataType       = 'erp'; opts.bands          = {''};
% opts.dataType       = 'power'; opts.bands          = {'hgam'};
% %opts.dataType       = 'power'; opts.bands          = {'delta','theta','alpha','beta','lgam','hgam'};
% opts.timeType       = 'Bin';
% opts.channelGroupingType      = 'channel';
% opts.timeFeatures   = 'trial';
% opts.extStr         = 'liblinearS0';
% opts.analysis       = 'trial';
% 
% % load decoding structure
% dataPath = ['../Results/Classification/group/' opts.dataType ...
%     '/' opts.channelGroupingType '/'];
% 
% fileName = ['SumallSubjsClass' opts.lockType 'Lock' opts.dataType cell2mat(opts.bands) '_tF' opts.timeFeatures '_tT' ...
%     opts.timeType '_gT' opts.channelGroupingType '_Solver' opts.extStr];
% 
% opts.fileName = fileName;
% load([dataPath fileName])
% 
% % load data structure
% dataPath = '../Results/';
% if strcmp(opts.dataType,'erp')
%     dataPath = [dataPath 'ERP_Data/group/'];
%     fileName = [opts.hems 'ERPsGroup' opts.lockType 'LocksubAmp' ...
%         opts.reference num2str(opts.nRefChans)];
% elseif strcmp(opts.dataType,'power')
%     dataPath = [dataPath 'Spectral_Data/group/'];
%     fileName = ['allERSPs' cell2mat(opts.bands) 'Group' opts.lockType  ...
%         'LocksublogPower' opts.reference num2str(opts.nRefChans)];
% end
% load([dataPath fileName])
% 
% out = matchChannelPatterns(data,S,opts);
% 
% savePath = '../Results/lagAnalysis/';
% if ~exist(savePath,'dir'), mkdir(savePath), end;
% fileName = ['lag' opts.analysis 'Analysis' opts.lockType 'Lock' opts.dataType cell2mat(opts.bands)];
% save([savePath fileName])
