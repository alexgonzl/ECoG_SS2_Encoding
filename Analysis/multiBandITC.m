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
trialSamps      = epochTime>=dur(1) & epochTime<=dur(2);
nTrialSamps     =sum(trialSamps);
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
for ch = 1:nChans
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
out.nTrialSamps = nTrialSamps;
out.trPh        = trPh;
out.nFreqs      = nFreqs;
out.behavior    = data.behavorInfo;
out.rois        = data.rois;

% behavior and conditions
out.conds{1}         = out.behavior.cond==1  ;  % objective abstract
out.conds{2}         = out.behavior.cond==2  ;  % objective concrete
out.conds{3}         = out.behavior.sResp==1 ;  % abstract response
out.conds{4}         = out.behavior.sResp==2 ;  % concrete response
out.conds{5}         = out.behavior.tResp<=2 ;  % old response at test
out.conds{6}         = out.behavior.tResp>=2 ;  % new response at test
out.conds{7}         = out.behavior.subRem==1; % remembered
out.conds{8}         = out.behavior.subForg==1; % forgotten
out.conds{9}         = (out.conds{1} & out.conds{3}) | ...
    (out.conds{2} & out.conds{4});              % correct semantic desc
out.conds{10}        = out.conds{9} & out.conds{7}; % correct at study and test
out.conds{11}        = out.behavior.studyRTs ;
out.conds{12}        = out.behavior.testRTs ;
out.conds{13}        = out.conds{7}&~isnan(out.conds{11}); % valid encoding, remembered trials

out.studyIDatTest   = out.behavior.studyIDatTest;
out.testRTs         = out.behavior.testRTs;
out.studyRTs        = data.behavior.studyRTs;

out.condOfInterest = 13; %#13: valid encoding. trial that were remembered.

% specific trial analyses
trials = out.conds{out.condOfInterest} ;
out.nTrials = sum(trials);
out.ITC     = abs(squeeze(mean(exp(1j*trPh(:,:,trials,:)),3)));
out.ITP     = angle(squeeze(mean(exp(1j*trPh(:,:,trials,:)),3)));

out.studyRTsToPhaseCorr = zeros(nChans,nFreqs,nTrialSamps);
out.testRTsToPhaseCorr  = zeros(nChans,nFreqs,nTrialSamps);

out.studyRTsPhaseGLMsPCA_R2 = zeros(nChans,nFreqs);
out.testRTsPhaseGLMsPCA_R2  = zeros(nChans,nFreqs);
out.studyRTsPhaseR2         = zeros(nChans,1);
out.testRTsPhaseR2          = zeros(nChans,1);
rts1   = out.studyRTs(trials);
rts2   = out.testRTs(trials);
% PCA-GLM approach.
for ch = 1:nChans
    for ff = 1:nFreqs
        alpha = squeeze(trPh(ch,ff,trials,:));
        out.studyRTsToPhaseCorr(ch,ff,:) = ccorr(rts1,alpha);
        out.testRTsToPhaseCorr(ch,ff,:)  = ccorr(rts2,alpha);
        
        v1=pca(cos(alpha)','NumComponents',12);
        m=fitglm(v1,rts1);
        out.studyRTsPhaseGLMsPCA_R2(ch,ff) = m.Rsquared.Ordinary;
        m=fitglm(v1,rts2);
        out.testRTsPhaseGLMsPCA_R2(ch,ff) = m.Rsquared.Ordinary;
    end
    alpha = squeeze(trPh(ch,:,trials,:));
    alpha2=permute(alpha,[2 1 3]); alpha2=alpha2(:,:);
    v1=pca(cos(alpha2)','NumComponents',12);
    m=fitglm(v1,rts1);
    out.studyRTsPhaseR2(ch) = m.Rsquared.Ordinary;
    m=fitglm(v1,rts2);
    out.testRTsPhaseR2(ch) = m.Rsquared.Ordinary;
    fprintf(' Analyses for chan %i completed',ch)
end
end
%% aux function
function r=ccorr(x,alpha)
nTr = size(x,1); % linear variable
nSp = size(alpha,2); % circular variable

r = zeros(nSp,1);
for sa = 1:nSp
    r(sa) = circ_corrcl(alpha(:,sa),x);
end
end


