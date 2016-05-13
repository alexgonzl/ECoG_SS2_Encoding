function out = multiBandITC(data,opts)
% function filters signal in single hertz bands epochs and computes the ITC
% across LPC channels

switch  data.lockType
    case 'preStim'
        data.trialDur = [-1 1]; dur = data.trialDur;
    case 'preStim2'
        data.trialDur = [-1.5 1.5]; dur = data.trialDur;
    case 'stim'
        data.trialDur = [-0.2 1.5]; dur = data.trialDur;
    case 'RT'
        data.trialDur = [-1.2 0.5]; dur = data.trialDur;
end

SR  = data.SR;
if isfield(opts,'downsample')
    downSampRate = opts.downsample;
else
    downSampRate = 4;
end
SR = SR/downSampRate;

data.nTrials        = numel(data.trialOnsets);
epochDur = [-2 2]; % before converting into trials and baseline correcting

% add behavioral data
data.behavior = data.behavorInfo; % 

% store the corresponding roi channel ids in data struct.
data.rois       = data.chanInfo;

% time is the total time around the trial that was stored in the data
% this is might be greater than the trial duration of interest
epochTime = linspace(epochDur(1),epochDur(2),ceil(diff(epochDur)*SR)); nEpSamps = numel(epochTime);

% determine the samples that correspond to the specified trial duration
trialSamps     = epochTime>=dur(1) & epochTime<=dur(2);

% Determine time samples in the trial
trialTime      = epochTime(trialSamps); data.trialTime = trialTime;

evOnsets = data.trialOnsets;
nEvents = numel(evOnsets);
epSamps = floor(epochTime*SR);
evIdx = zeros(nEvents,nEpSamps);

switch  data.lockType    
    case 'RT'
        offset = floor(data.behavior.studyRTs*SR);
    otherwise
        offset = zeros(nEvents,1);
end

% epoch samples matrix
for ev = 1:nEvents
    evIdx(ev,:) = epSamps+evOnsets(ev)+offset(ev);
end

% LPC channels
channels = data.rois.LPCchannels;
nChannels = numel(channels);

%
validTrials = find(~isnan(offset))';
X = data.amp(channels,:);data.signal=[]; data.amp=[]; data.phase=[];
Y = multiBandFilter(X,maxF,SR)

% low pass the amplitude
%X=channelFilt(X,SR,20,[],[]);


erp = nan(data.nChans,nEvents,nEpSamps);
for ch = 1:data.nChans
    x = X(ch,:);
    for tr = validTrials
        erp(ch,tr,:) = x(evIdx(tr,:));
    end
end

% Percent signal change for amplitude.
switch data.analysisType
    case 'Power'
        erp = erp.^2;
    case 'logPower'
        erp = 20*log10(abs(erp));
end

% Correct for baseline fluctuations before trial onset. In case of RT
% locked analysis, it uses the baselines from the stim locked (loaded outside)
if ~strcmp(data.lockType,'RT')
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

end




