function data=calcERP(data)
% function epochs the raw signals, separates by trial condition and
% calculates simple univariate stats
%
% data should be the output of reReferenceData.m
%
% data dependecies:
%   loads behavioral data
%
% file dependencies:
%   channelFilt.m
%   CalculateBadTrials.m


data.lowpass = 20; lp = data.lowpass; SR = data.SR;
switch  data.lockType
    case 'preStim'
        data.trialDur = [-1 1]; dur = data.trialDur;
        data.baseLine = [-1 1];
    case 'stim'
        data.trialDur = [-0.2 1.5]; dur = data.trialDur;
        data.baseLine = [-0.2 0];
    case 'RT'
        data.trialDur = [-1 0.2]; dur = data.trialDur;
        % baseline period used from the stim locked data 
end

data.nTrials        = numel(data.trialOnsets);
epochDur = [-2 2]; % before converting into trials and baseline correcting

% add behavioral data
data.behavior = data.behavorInfo; % hack... fix 5/7/14

% store the corresponding roi channel ids in data struct.
data.rois       = data.chanInfo;

% time is the total time around the trial that was stored in the data
% this is might be greater than the trial duration of interest
epochTime = linspace(epochDur(1),epochDur(2),ceil(diff(epochDur)*data.SR)); nEpSamps = numel(epochTime);

% determine the samples that correspond to the specified trial duration
trialSamps     = epochTime>=dur(1) & epochTime<=dur(2);

% Determine time samples in the trial
trialTime      = epochTime(trialSamps); data.trialTime = trialTime;

evOnsets = data.trialOnsets; % in samples
nEvents = numel(evOnsets);
epSamps = floor(epochTime*data.SR);
evIdx = zeros(nEvents,nEpSamps);

switch  data.lockType
    case {'stim','preStim'}
        offset = zeros(nEvents,1);
    case 'RT'
        offset = floor(data.behavior.allRTs*data.SR);
end

for ev = 1:nEvents
    evIdx(ev,:) = epSamps+evOnsets(ev)+offset(ev);
end

validTrials = find(~isnan(offset))';

% low pass the data
X = data.signal;data.signal=[];
X = channelFilt(X,SR,lp,[],[]);

% epoch the data
erp = nan(data.nChans,nEvents,nEpSamps);
for ch = 1:data.nChans
    x = X(ch,:);
    for tr = validTrials
        erp(ch,tr,:) = x(evIdx(tr,:));
    end
end

data.evIdx = evIdx;

% Percent signal change for amplitude.
switch data.analysisType
    case 'Power'
        erp = erp.^2;
    case 'logPower'
        erp = 20*log10(erp);
    otherwise 
        % don't do anything
end

% Correct for baseline fluectuations before trial onset. In case of RT
% locked analysis, it uses the baselines from the stim locked (loaded outside)
if  ~strcmp(data.lockType,'RT')
    baselineIdx = epochTime<=data.baseLine(2) & epochTime>= data.baseLine(1);
        data.baseLineMeans = nanmean(erp(:,:,baselineIdx),3);
end

switch data.baselineType
    case 'sub'
        erp = bsxfun(@minus,erp,data.baseLineMeans);
    case 'rel'
        erp = bsxfun(@rdivide,erp,data.baseLineMeans);
end

% Take relevant samples of trial
erp = erp(:,:,trialSamps); data.erp = erp;

% Compute the trial evoked mean and variance by channel
evokedSamples=data.trialTime>0;
data.evokedMean = zeros(data.nChans, data.nTrials);
data.evokedVar = zeros(data.nChans, data.nTrials);
for ch = 1:data.nChans
    data.evokedMean(ch,:) = mean(data.erp(ch,:,evokedSamples),3);
    data.evokedVar(ch,:)  = var(data.erp(ch,:,evokedSamples),0,3);
end

end




