function data=calcERSP(data)
% function epochs the raw signals, separates by trial condition and
% calculates simple univariate stats

switch  data.lockType
    case 'preStim'
        data.trialDur = [-1 1]; dur = data.trialDur;
        data.baseLine = [-1 1];
    case 'preStim2'
        data.trialDur = [-1.5 1.5]; dur = data.trialDur;
        data.baseLine = [-1.5 -1];
    case 'stim'
        data.trialDur = [-0.2 1.5]; dur = data.trialDur;
        data.baseLine = [-0.2 0];
    case 'RT'
        data.trialDur = [-1.2 0.5]; dur = data.trialDur;
        % baseline period used from the stim locked data
end

data.nTrials        = numel(data.trialOnsets);
epochDur = [-2 2]; % before converting into trials and baseline correcting

% add behavioral data
data.behavior = data.behavorInfo; % 

% store the corresponding roi channel ids in data struct.
data.rois       = data.chanInfo;

% time is the total time around the trial that was stored in the data
% this is might be greater than the trial duration of interest
epochTime = linspace(epochDur(1),epochDur(2),ceil(diff(epochDur)*data.SR)); nEpSamps = numel(epochTime);

% determine the samples that correspond to the specified trial duration
trialSamps     = epochTime>=dur(1) & epochTime<=dur(2);

% Determine time samples in the trial
trialTime      = epochTime(trialSamps); data.trialTime = trialTime;

evOnsets = data.trialOnsets;
nEvents = numel(evOnsets);
epSamps = floor(epochTime*data.SR);
evIdx = zeros(nEvents,nEpSamps);

switch  data.lockType    
    case 'RT'
        offset = floor(data.behavior.studyRTs*data.SR);
    otherwise
        offset = zeros(nEvents,1);
end

for ev = 1:nEvents
    evIdx(ev,:) = epSamps+evOnsets(ev)+offset(ev);
end

validTrials = find(~isnan(offset))';
X = data.amp;data.signal=[]; data.amp=[]; data.phase=[];
% low pass the amplitude
X=channelFilt(X,data.SR,20,[],[]);


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




