function data = calcMI(AmpData,PhaseData)
% dependencies:
%   subjChanInfo
%   subjExptInfo
%   modIndex

data = [];
data.ampSR              = AmpData.SR;
data.phSR               = PhaseData.SR;
data.rois               = PhaseData.rois;
data.chans              = data.rois.LPC;
data.nChans             = numel(data.chans);
data.nFreqs             = PhaseData.nFreqs;
data.downRate           = data.ampSR/data.phSR;
data.behavior           = AmpData.behavior;
data.conds              = PhaseData.conds;
data.nTrials            = AmpData.nTrials;
data.nPhTrSamps         = PhaseData.nTrialSamps;

if any(strcmp(AmpData.lockType,{'stim','preStim2'}))
    ampValidSamples   = AmpData.trialTime>AmpData.baseLine(2);   
    phaseValidSamples = PhaseData.trialTime>AmpData.baseLine(2);          
else
    ampValidSamples   = true(size(AmpData.trialTime));
    phaseValidSamples = true(size(PhaseData.trialTime));
end

% get amplitude
amp = AmpData.erp(data.chans,:,ampValidSamples);

% Separate trials based on condition
data.nPhBins = 20;      
data.phaseBinEdges      = linspace(-pi,pi,data.nPhBins+1); 
phaseBinEdges           = data.phaseBinEdges;
data.phaseBinCenters    = mean([phaseBinEdges(1:end-1)' phaseBinEdges(2:end)'],2);

data.MI = zeros(data.nChans,data.nFreqs,data.nTrials,data.nPhBins);
data.mMagVec = zeros(data.nChans,data.nFreqs,data.nTrials);
data.mPhVec = zeros(data.nChans,data.nFreqs,data.nTrials);
for cc = 1:data.nChans
    % amplitude
    temp   = squeeze(amp(cc,:,:));
    x = downsample(temp',data.downRate)';

    % throw away end samples to match dimensions.
    dsamp = size(x,2)-sum(phaseValidSamples);
    if dsamp>0 && dsamp<data.downRate
        x(:,dsamp) = [];
    elseif dsamp<0
        error('more phase samples than amplitude samples');
    elseif dsamp>=data.downRate
        error('matching phase and amplitude not possible with current sampling')        
    end

    for ff = 1:data.nFreqs
        phase   = squeeze(PhaseData.trPh(cc,ff,:,phaseValidSamples));
        cMI     = zeros(data.nTrials,data.nPhBins);
        for tr = 1:data.nTrials
            cMI(tr,:) = modIndex(phase(tr,:),x(tr,:),phaseBinEdges);
        end
        data.MI(cc,ff,:,:) = cMI;
        
        yy=mean(x.*exp(1j*phase),2);
        data.mMagVec(cc,ff,:) =  abs(yy);
        data.mPhVec(cc,ff,:)  =  angle(yy);
    end  
end

return