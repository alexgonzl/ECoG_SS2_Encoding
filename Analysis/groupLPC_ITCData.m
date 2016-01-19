function data = groupLPC_ITCData(opts)
% groups LPC channel data by subject.

% data parameters
subjects    = opts.subjects;
reference   = opts.reference;
lockType    = opts.lockType;
band        = opts.band;

type        = ['/ITC_Data/' band '/'];
preFix      = ['ITC' band];
extension   = [lockType 'sublogPower' reference] ;
fileName    = [preFix extension];

nSubjs  = numel(subjects);

data                = [];
data.options        = opts;
data.prefix         = preFix;
data.extension      = extension;
data.options        = opts;
data.subjChans      = [];
data.LPCchanId      = [];
data.hemChan        = [];
data.ROIid          = []; data.ROIs = {'IPS','SPL','AG'};
data.subROIid       = []; data.subROIs = {'pIPS','aIPS','pSPL','aSPL'};

data.nRTsplits      = 3;
data.ITC_testRTSplits = [];
data.ITC_studyRTSplits= [];

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
    
    chans = temp.LPC;
    data.phase{s}           = squeeze(dataIn.data.phaseResp(chans,:,:));
    
    data.ITC_testRTSplits = cat(2,data.ITC_testRTSplits,squeeze(dataIn.data.ITC_testRTSplits(:,chans,:)));
    data.ITC_studyRTSplits =cat(2,data.ITC_studyRTSplits,squeeze(dataIn.data.ITC_studyRTSplits(:,chans,:)));
    
    data.testRTquantiles{s} = dataIn.data.testRTquantiles;
    data.studyRTquantiles{s}= dataIn.data.studyRTquantiles;    
    data.studyIDatTest{s}   = dataIn.data.behavior.studyIDatTest;
    data.testRTs{s}         = dataIn.data.behavior.testRTs;
    data.studyRTs{s}        = dataIn.data.behavior.studyRTs;
    
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

data.BinITC_testRTSplits     = zeros(data.nRTsplits,data.nChans,data.nBins);
data.BinITC_studyRTSplits    = zeros(data.nRTsplits,data.nChans,data.nBins);


for ch = 1:data.nChans
    for sp = 1:data.nRTsplits
        data.BinITC_testRTSplits(sp,ch,:) = binTrials(data.ITC_testRTSplits(sp,ch,:),data.BinSamps);
        data.BinITC_studyRTSplits(sp,ch,:)= binTrials(data.ITC_studyRTSplits(sp,ch,:),data.BinSamps);
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

