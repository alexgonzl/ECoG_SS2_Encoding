function data = groupLPCData(opts)
% groups LPC channel data by subject.
% takes outputs from calcERP or calcERSP
%
% dependencies:
%       getBinSamps
%       signtest2
%       ranksum2

% data parameters
subjects    = opts.subjects;
reference   = opts.reference;
lockType    = opts.lockType;
band        = opts.band;

% erps or a spectral band
switch band
    case 'erp'
        baselineType = 'sub';
        analysisType = 'Amp';
        type    = 'ERP_Data/';
        preFix = 'ERPs' ;
    otherwise
        baselineType = 'sub';
        analysisType = 'logPower';
        type    = ['Spectral_Data/' band '/'];
        preFix = ['ERSPs' band];
end

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
data.ROIid          = []; data.ROIs = {'IPS','SPL','AG'};
data.subROIid       = []; data.subROIs = {'pIPS','aIPS','pSPL','aSPL'};
data.RTquantiles    = [0.1:0.1:0.9]; % quanntiles for analysis

data.condNames      = {'abstract','concrete','absResp','conResp',...
    'oldResp','newResp', 'Rem','Forg', ...
    'studyRTs','testRTs'};

for s = 1:nSubjs
    
    dataIn = load([opts.dataPath subjects{s} '/' type fileName '.mat']);
    
    temp                    = subjChanInfo(subjects{s});
    data.subjHemsIds{s}        = temp.hemisphere;
    data.nSubjLPCchans(s)   = numel(temp.LPC);
    data.LPCchanId          = [data.LPCchanId; temp.LPC'];
    data.subjChans          = [data.subjChans; s*ones(data.nSubjLPCchans(s),1)];
    data.hemChan            = [data.hemChan; (strcmp(data.subjHemsIds{s},'r')+1)*ones(data.nSubjLPCchans(s),1)]; % one for lefts
    
    ROIid = 1*ismember(temp.LPC,temp.IPS) + ...
        2*ismember(temp.LPC,temp.SPL) + ...
        3*ismember(temp.LPC,temp.AG);
    subROIid = 1*ismember(temp.LPC,temp.pIPS) + ...
        2*ismember(temp.LPC,temp.aIPS) + ...
        3*ismember(temp.LPC,temp.pSPL) + ...
        4*ismember(temp.LPC,temp.aSPL);
    
    data.ROIid              = [data.ROIid; ROIid'];
    data.subROIid           = [data.subROIid; subROIid'];
    
    data.ERP{s}             = squeeze(dataIn.data.erp(temp.LPC,:,:));
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
data.winSize        = 0.10; % in seconds
data.sldWin         = data.winSize;
[data.BinSamps , data.Bins] = getBinSamps(data.winSize,data.sldWin,data.trialTime);
data.nBins          = size(data.Bins,1);

data.nChans                 = sum(data.nSubjLPCchans);

data.BinERP                 = cell(nSubjs,1);
data.meanChResp             = zeros(data.nChans,data.nTrialSamps);  
data.meanBinChResp          = zeros(data.nChans,data.nBins);  
data.nTrialsInAnalysis      = zeros(nSubjs,1);
data.dataToStudyRTsCorr     = zeros(data.nChans,data.nBins);
data.dataToStudyRTsCorrP    = zeros(data.nChans,data.nBins);
data.dataToTestRTsCorr      = zeros(data.nChans,data.nBins);
data.dataToTestRTsCorrP     = zeros(data.nChans,data.nBins);


ch = 1;
for s = 1:nSubjs
    % remembered items and correctly identified semantics
    trials = data.conds{s,10};
    nChans = data.nSubjLPCchans(s);
    nTrials = numel(trials);
    %data.nTrialsInAnalysis(s) = nTrials;
    
    testRTs = data.testRTs{s};
    studyRts = data.studyRTs{s};
    
    % pre-allocate
    data.BinERP{s}             = zeros(nChans,nTrials,data.nBins);
    
    for Sch = 1: data.nSubjLPCchans(s)
        Z       = squeeze(data.ERP{s}(Sch,:,:));
        data.BinERP{s}(Sch,:,:) = binTrials(Z,data.BinSamps);
        
        data.meanChResp(ch,:)   = mean(Z);
        data.meanBinChResp(ch,:) =mean(data.BinERP{s}(Sch,:,:));
        
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
