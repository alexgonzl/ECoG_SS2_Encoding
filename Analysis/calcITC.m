function data = calcITC(data)
% data = calcITC(data)
% function epochs the phase signals, separates by trial condition computes
% inter trial coherence and permutation tests to quantify significance per
% channel
%
% dependencies:
%   subjChanInfo
%   subjExptInfo
%   itc

data.nPhBins = 200;
data.nBoot = 20;
data.nShuf = 100;
data.RTsplits = 3;

switch  data.lockType
    case 'stim'
        data.trialDur = [-0.2 1.5]; dur = data.trialDur;
        data.baseLine = [-0.2 0];
    case 'preStim'
        data.trialDur = [-1 1]; dur = data.trialDur;
        data.baseLine = [-1 1];
    case 'RT'
        data.trialDur = [-1 0.2]; dur = data.trialDur;
        % baseline period used from the stim locked data
end

data.nTrials        = numel(data.trialOnsets);
epochDur = data.trialDur; % before converting into trials and baseline correcting

% add behavioral data
data.behavior = data.behavorInfo; 

data.rois = data.chanInfo;

% time is the total time around the trial that was stored in the data
% this is might be greater than the trial duration of interest
epochTime = linspace(epochDur(1),epochDur(2),round(diff(epochDur)*data.SR)); nEpSamps = numel(epochTime);

% determine the samples that correspond to the specified trial duration
trialSamps     = epochTime>=dur(1) & epochTime<=dur(2); 

% Determine time samples in the trial
trialTime      = epochTime(trialSamps); data.trialTime = trialTime;

% get the reation times from behavioral data
data.studyRTs = data.behavior.studyRTs;
data.studyRTquantiles = quantile(data.studyRTs,(1:data.RTsplits)/data.RTsplits);
data.testRTs  = data.behavior.testRTs;
data.testRTquantiles = quantile(data.testRTs,(1:data.RTsplits)/data.RTsplits);

evOnsets = data.trialOnsets;
nEvents = numel(evOnsets);
epSamps = floor(epochTime*data.SR);
evIdx = zeros(nEvents,nEpSamps);

switch  data.lockType
    case {'stim','preStim'}
        offset = zeros(nEvents,1);
    case 'RT'
        offset = round(data.behavior.studyRTs*data.SR);
end

for ev = 1:nEvents
    evIdx(ev,:) = epSamps+evOnsets(ev)+offset(ev);
end

validTrials = find(~isnan(offset))';

    
% trials to be analyzed.
data.correctSemanticResp = (data.behavior.sResp==1 & data.behavior.tResp==1) |...
    (data.behavior.sResp==2 & data.behavior.tResp==2) ;
data.correctReps   = data.correctSemanticResp & data.behavior.subRem;
data.correctRepsIDs = find(data.correctReps);

% store the phase;
nChans  = data.nChans; nShuf = data.nShuf; nBoot = data.nBoot;
X = data.phase; data.signal=[]; data.amp=[]; data.phase=[];
phase      = nan(nChans,nEvents,nEpSamps);
for ch = 1:nChans
    x = X(ch,:);
    for tr = validTrials
        phase(ch,tr,:) = x(evIdx(tr,:));
    end
end
data.phaseResp = phase; 

% get the reation times from behavioral data
data.studyRTs = data.behavior.studyRTs;
data.studyRTquantiles = quantile(data.studyRTs(data.correctReps),(1:data.RTsplits)/data.RTsplits);
data.testRTs  = data.behavior.testRTs;
data.testRTquantiles = quantile(data.testRTs(data.correctReps),(1:data.RTsplits)/data.RTsplits);

data.testTrialSets = cell(data.RTsplits,1);
data.studyTrialSets = cell(data.RTsplits,1);
for sp = 1:data.RTsplits
    if sp==1
        temp1 =  data.testRTs<=data.testRTquantiles(sp);
        temp2 =  data.studyRTs<=data.studyRTquantiles(sp);
    else 
        temp1 =  data.testRTs>data.testRTquantiles(sp-1) & data.testRTs<=data.testRTquantiles(sp);
        temp2 =  data.studyRTs>data.studyRTquantiles(sp-1) & data.studyRTs<=data.studyRTquantiles(sp);
    end
        data.testTrialSets{sp} = temp1 & data.correctReps;
        data.studyTrialSets{sp} = temp2 & data.correctReps;
end


data.ITC_testRTSplits  = zeros(data.RTsplits,data.nChans,nEpSamps);
data.ITC_BootSD_testRTSplits  = zeros(data.RTsplits,data.nChans,nEpSamps);
data.ITC_studyRTSplits = zeros(data.RTsplits,data.nChans,nEpSamps);
data.ITC_BootSD_studyRTSplits = zeros(data.RTsplits,data.nChans,nEpSamps);

% Get ITC values per RT splits
for ch = 1:nChans
    for sp = 1:data.RTsplits
        X = squeeze(data.phaseResp(ch,data.testTrialSets{sp},:));
        data.ITC_testRTSplits(sp,ch,:)= itc(X);
        data.ITC_BootSD_testRTSplits(sp,ch,:) = std(bootstrp(nBoot,@itc,X));
        
        X = squeeze(data.phaseResp(ch,data.studyTrialSets{sp},:));
        data.ITC_studyRTSplits(sp,ch,:)= itc(X);
        data.ITC_BootSD_studyRTSplits(sp,ch,:) = std(bootstrp(nBoot,@itc,X));
    end
end

