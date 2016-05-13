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
if isfield(opts,'nFreqs')
    nFreqs = opts.nFreqs;
else
    nFreqs = 30;
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
trialTime      = epochTime(trialSamps); 

evOnsets = data.trialOnsets; % in samples at downsampled rate
evOnsets = round(evOnsets/downSampRate);
nEvents = numel(evOnsets);
epSamps = floor(epochTime*SR); % epoch time in samples 
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
validTrials = find(~isnan(offset))';

% LPC channels
chans   = data.rois.LPC;
nChans  = numel(chans);

% decimate data
X = data.signal(chans,:);
x = decimate(X(1,:),downSampRate);
nSamps  = numel(x);
Y       = zeros(nChans,nSamps);
Y(1,:)  = x;
for ch = 2:nChans
    Y(ch,:) = decimate(X(ch,:),downSampRate);
end
clear X;

Z = multiBandFilter(Y,nFreqs,SR);
clear Y;

% decompose and extract phase per channel and frequency
Ang = zeros(size(Z));
parfor ch = 1:nChans    
    [~,Ang(ch,:,:)] = multiBandDecomp(squeeze(Z(ch,:,:)));
end
clear Z;

trPh = nan(nChans,nFreqs,nEvents,nEpSamps);
for ch = 1:nChans
    for ff = 1:nFreqs
        x = squeeze(Ang(ch,ff,:,:));
        for tr = validTrials
            trPh(ch,ff,tr,:) = x(evIdx(tr,:));
        end
    end
end

% Take relevant samples of trial
trPh = trPh(:,:,:,trialSamps); 


%% output
out             = [];
out.SR          = SR;
out.nEvents     = nEvents;
out.trialDur    = dur;
out.trialTime   = trialTime;
out.trPh        = trPh;
out.nFreqs      = nFreqs;
out.behavior    = data.behavorInfo;
out.rois        = data.rois;



