function data = groupMultiBandITCData(opts)
% groups LPC channel data by subject.
% takes outputs from multiBandITC
%
% dependencies:
%       getBinSamps
%       signtest2
%       ranksum2

% data parameters
nFreq       = opts.nFreq;
subjects    = opts.subjects;
reference   = opts.reference;
lockType    = opts.lockType;

extension = [lockType baselineType analysisType reference] ;
fileName = [preFix extension];

nSubjs  = numel(subjects);

data                = [];
data.options        = opts;
data.prefix         = preFix;
data.extension      = extension;
data.options        = opts;
data.subjChans      = [];
data.ERP            = [];
data.LPCchanId      = [];
data.hemChan        = [];
data.ROIid          = []; data.ROIs = {'IPS','SPL','AG','TPJ','SMG'};
data.subROIid       = []; data.subROIs = {'pIPS','aIPS','pSPL','aSPL'};
data.RTquantiles    = 0.1:0.1:0.9; % quanntiles for analysis

data.condNames      = {'abstract','concrete','absResp','conResp',...
    'oldResp','newResp', 'Rem','Forg', 'correctSemantic', 'correctStudyTest',...
    'studyRTs','testRTs','validEncRemem'};

for s = 1:nSubjs
    
    % load subjects data
    dataIn = load([opts.dataPath subjects{s} '/' type fileName '.mat']);
    
    % select, ID and store relevant ROI channels
    temp                    = subjChanInfo(subjects{s});
    data.subjHemsIds{s}        = temp.hemisphere;
    data.nSubjLPCchans(s)   = numel(temp.LPC);
    data.LPCchanId          = [data.LPCchanId; temp.LPC'];
    data.subjChans          = [data.subjChans; s*ones(data.nSubjLPCchans(s),1)];
    data.hemChan            = [data.hemChan; (strcmp(data.subjHemsIds{s},'r')+1)*ones(data.nSubjLPCchans(s),1)]; % one for lefts
    
    ROIid = 1*ismember(temp.LPC,temp.IPS) + ...
        2*ismember(temp.LPC,temp.SPL) + ...
        3*ismember(temp.LPC,temp.AG) + ... 
        4*ismember(temp.LPC,temp.SMG) + ...
        5*ismember(temp.LPC,temp.TPJ) ;
    
    subROIid = 1*ismember(temp.LPC,temp.pIPS) + ...
        2*ismember(temp.LPC,temp.aIPS) + ...
        3*ismember(temp.LPC,temp.pSPL) + ...
        4*ismember(temp.LPC,temp.aSPL);
    
    data.ROIid              = [data.ROIid; ROIid'];
    data.subROIid           = [data.subROIid; subROIid'];
    
    % get data for LPC channels
    data.ERP{s}             = squeeze(dataIn.data.erp(temp.LPC,:,:));
    
    % behavior and conditions
    data.behavior{s}        = dataIn.data.behavior;    
    data.conds{s,1}         = data.behavior{s}.cond==1  ;  % objective abstract
    data.conds{s,2}         = data.behavior{s}.cond==2  ;  % objective concrete
    data.conds{s,3}         = data.behavior{s}.sResp==1 ;  % abstract response
    data.conds{s,4}         = data.behavior{s}.sResp==2 ;  % concrete response
    data.conds{s,5}         = data.behavior{s}.tResp<=2 ;  % old response at test
    data.conds{s,6}         = data.behavior{s}.tResp>=2 ;  % new response at test
    data.conds{s,7}         = data.behavior{s}.subRem==1; % remembered
    data.conds{s,8}         = data.behavior{s}.subForg==1; % forgotten
    data.conds{s,9}         = (data.conds{s,1} & data.conds{s,3}) | ...
        (data.conds{s,2} & data.conds{s,4});              % correct semantic desc
    data.conds{s,10}        = data.conds{s,9} & data.conds{s,7}; % correct at study and test
    data.conds{s,11}        = data.behavior{s}.studyRTs ;
    data.conds{s,12}        = data.behavior{s}.testRTs ;
    data.conds{s,13}        = data.conds{s,7}&~isnan(data.conds{s,11}); % valid encoding, remembered trials
    
    data.studyIDatTest{s}   = data.behavior{s}.studyIDatTest;
    data.testRTs{s}         = dataIn.data.behavior.testRTs;
    data.studyRTs{s}        = dataIn.data.behavior.studyRTs;
    data.testRTsQuants{s}   = quantile(data.testRTs{s}(data.conds{s,7}==1),data.RTquantiles);
    
end

% trial time info
data.baselineType   = dataIn.data.baselineType;
data.SR             = dataIn.data.SR;
data.nTrialSamps    =numel(dataIn.data.trialTime);
data.trialDur       = dataIn.data.trialDur;
data.trialTime      = linspace(data.trialDur(1),data.trialDur(2),data.nTrialSamps);
% bin info
data.winSize                = 0.05; % in seconds
data.sldWin                 = data.winSize;
% all bins
[data.BinSamps , data.Bins] = getBinSamps(data.winSize,data.sldWin,data.trialTime);
data.nBins                  = size(data.Bins,1);

% Get samples to be used in analyses. Excludes baseline period.
% RT-locked data has baseline outside the epochs. preStim has no baseline.
% other cases should include the samples greater than the last sample
% included in the baseline
if any(strcmp(opts.lockType,{'RT','preStim'}))
    data.validAnalysisSamples = true(size(data.trialTime));
    % analysis bins
    data.AnalysisBins = true(size(data.Bins(:,2),1),1);
else
    data.baseLine       = dataIn.data.baseLine;
    data.validAnalysisSamples = data.trialTime>data.baseLine(2);   
    data.AnalysisBins         = data.Bins(:,2)>=data.baseLine(2);
end
data.nValidAnalysisSamps = sum(data.validAnalysisSamples);

% number of channels
data.nChans                 = sum(data.nSubjLPCchans);

% Pre-allocations
data.BinERP                 = cell(nSubjs,1);

% all these correspond to the condition of interest: #13: valid encoding
% trial that were remembered.
data.condOfInterest         = 13;
data.mChResp                = zeros(data.nChans,data.nTrialSamps);  
data.tChResp                = zeros(data.nChans,data.nTrialSamps);
data.chScore                = zeros(data.nChans,1); % for valid samples
data.mBinChResp             = zeros(data.nChans,data.nBins);  
data.tBinChResp             = zeros(data.nChans,data.nBins);  
data.chBinScore             = zeros(data.nChans,1);

data.nTrialsInAnalysis      = zeros(nSubjs,1);
data.trialsInAnalysis       = cell(nSubjs,1);
data.dataToStudyRTsCorr     = zeros(data.nChans,data.nBins);
data.dataToStudyRTsCorrP    = zeros(data.nChans,data.nBins);
data.dataToTestRTsCorr      = zeros(data.nChans,data.nBins);
data.dataToTestRTsCorrP     = zeros(data.nChans,data.nBins);


ch = 1;
for s = 1:nSubjs
    % remembered items and correctly identified semantics
    trials = data.conds{s,data.condOfInterest };
    nChans = data.nSubjLPCchans(s);
    nTrials = numel(trials);
    data.nTrialsInAnalysis(s) = sum(trials);
    data.trialsInAnalysis{s} = trials;
    
    testRTs = data.testRTs{s};
    studyRts = data.studyRTs{s};
    
    % pre-allocate
    data.BinERP{s}             = zeros(nChans,nTrials,data.nBins);
    
    for Sch = 1: data.nSubjLPCchans(s)
        Z       = squeeze(data.ERP{s}(Sch,:,:));
        data.BinERP{s}(Sch,:,:) = binTrials(Z,data.BinSamps);
        
        % summaries only for trials of interest
        % Continous data        
        data.mChResp(ch,:)   = mean(Z(trials,:));
        [~,~,~,t]=ttest(Z(trials,:));
        data.tChResp(ch,:) = t.tstat;
        data.chScore(ch)   = norm(t.tstat(data.validAnalysisSamples),inf);        
        
        % Bin data
        data.mBinChResp(ch,:) =mean(data.BinERP{s}(Sch,trials,:));
        [~,~,~,t]=ttest(data.BinERP{s}(Sch,trials,:));
        data.tBinChResp(ch,:) = t.tstat;
        data.chBinScore(ch)   = norm(squeeze(t.tstat(data.AnalysisBins)),inf); 
        
        % compute relationship of data to RTs
        [data.dataToStudyRTsCorr(ch,:),data.dataToStudyRTsCorrP(ch,:)] = ...
            corrDataToRTs(squeeze(data.BinERP{s}(Sch,trials,:)),studyRts(trials));
        
        [data.dataToTestRTsCorr(ch,:),data.dataToTestRTsCorrP(ch,:)] = ...
            corrDataToRTs(squeeze(data.BinERP{s}(Sch,trials,:)),testRTs(trials));
        
        ch = ch + 1;
    end
end
end

function binZ = binTrials(Z,samps)
nT = size(Z,1);
nB = size(samps,1);

binZ = nan(nT,nB);
for tt = 1:nT
    for bb = 1:nB
        binZ(tt,bb) = mean(Z(tt,samps(bb,:)));
    end
end
end
