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

extension = [lockType 'Lock' baselineType analysisType reference] ;
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
    
    temp = subjChanInfo(subjects{s});
    data.nSubjLPCchans(s)   = numel(temp.LPC);
    data.LPCchanId          = [data.LPCchanId; temp.LPC'];
    data.subjChans          = [data.subjChans; s*ones(data.nSubjLPCchans(s),1)];
    data.hemChan            = [data.hemChan; (strcmp(opts.hemId{s},'r')+1)*ones(data.nSubjLPCchans(s),1)]; % one for lefts
    
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
    data.conds{s,6}         = data.behavior{s}.tResp>=2 ; % new response at test
    data.conds{s,7}         = data.behavior{s}.subRem==1 ; % remembered
    data.conds{s,8}         = data.behavior{s}.subForg==1; % forgotten
    data.conds{s,9}         = data.behavior{s}.studyRTs ;
    data.conds{s,10}        = data.behavior{s}.testRTs ;

    data.studyIDatTest{s}   = data.behavior{s}.studyIDatTest;
    data.testRTs{s}         = dataIn.data.behavior.testRTs;
    data.studyRTs{s}        = dataIn.data.behavior.studyRTs;
    data.testRTsQuants{s}   = quantile(data.testRTs{s}(data.conds{s,5}==1),data.RTquantiles);
    
end

nChans = numel(data.subjChans);

% trial time info
data.baselineType   = dataIn.data.baselineType;
data.SR             = dataIn.data.SR;
data.nTrialSamps    =numel(dataIn.data.trialTime);
data.trialDur       = dataIn.data.trialDur;
data.trialTime      = linspace(data.trialDur(1),data.trialDur(2),data.nTrialSamps);

% bin info
data.winSize        = 0.1; % in seconds
data.sldWin         = data.winSize;
[data.BinSamps , data.Bins] = getBinSamps(data.winSize,data.sldWin,data.trialTime);
data.nBins          = size(data.Bins,1);

ch = 1;
data.BinERP     = cell(nSubjs,1);
for s = 1:nSubjs
    for Sch = 1: data.nSubjLPCchans(s)
        Z       = squeeze(data.ERP{s}(Sch,:,:));
        data.BinERP{s}(Sch,:,:) = binTrials(Z,data.BinSamps);
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

function data = getEffectScores(Z,RTs,C1,C2,data,ch,preFix)
% data = getEffectScores(Z,RTs,H,CRs,data,preFix)
% inputs:
% Z     -> data matrix with 1st dimension of Ntrials. 2nd dim of timepoints
% RTs   -> vector with Ntrials
% C1    -> cond1, logical vector Ntrials
% C2   -> cond2, logical vector Ntrials
% data  -> data structure to write output to.
% preFix-> preFix string for the outputs

nC1    = sum(C1);
nC2    = sum(C2);

% original sampled data
X       = Z(C1,:);
Y       = Z(C2,:);

% store means for each condition
data.([preFix 'mCond1'])(ch,:)        = mean(X);
data.([preFix 'mCond2'])(ch,:)        = mean(Y);

if ~strcmp(preFix,'')
    % test each condition against 0
    [~,data.([preFix 'zCond1'])(ch,:)]                = signtest2(X);
    [~,data.([preFix 'zCond2'])(ch,:)]                = signtest2(Y);
    % test conditions against each other
    [data.([preFix 'PValZ'])(ch,:), data.([preFix 'ZStat'])(ch,:)] = ranksum2(X,Y);
    
    % correlate each condition to RT (transformed to z values)
    data.([preFix 'cCond1'])(ch,:)        = corr(X,log10(RTs(C1)));
    data.([preFix 'cCond2'])(ch,:)        = corr(Y,log10(RTs(C2)));
    data.([preFix 'ZcCond1'])(ch,:)       = atanh(data.([preFix 'cCond1'])(ch,:))*sqrt(nC1-3);
    data.([preFix 'ZcCond2'])(ch,:)       = atanh(data.([preFix 'cCond2'])(ch,:))*sqrt(nC2-3);
    data.([preFix 'ZcStat'])(ch,:)        = data.([preFix 'ZcCond1'])(ch,:)-data.([preFix 'ZcCond2'])(ch,:);
end
end