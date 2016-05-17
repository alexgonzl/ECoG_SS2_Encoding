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
data.RTquant            = [0.4 0.6];
data.studyRTs           = PhaseData.studyRTs;
data.testRTs            = PhaseData.testRTs;
data.studyRTsQuants     = PhaseData.studyRTsQuants;
data.testRTsQuants      = PhaseData.testRTsQuants;

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
data.nPhBins = 10;      
data.phaseBinEdges      = linspace(-pi,pi,data.nPhBins+1); 
phaseBinEdges           = data.phaseBinEdges;
data.phaseBinCenters    = mean([phaseBinEdges(1:end-1)' phaseBinEdges(2:end)'],2);

data.MI    = zeros(data.nChans,data.nFreqs,data.nTrials,data.nPhBins);
data.studyRTsMI_R2 = zeros(data.nChans,data.nFreqs);
data.studyRTsMI_R2_allF = zeros(data.nChans,1);
data.testRTsMI_R2 = zeros(data.nChans,data.nFreqs);
data.testRTsMI_R2_allF = zeros(data.nChans,1);

data.mMagVec = zeros(data.nChans,data.nFreqs,data.nTrials);
data.mPhVec = zeros(data.nChans,data.nFreqs,data.nTrials);
data.studyRTsmMagVec_R2_allF = zeros(data.nChans,1);
data.testRTsmMagVec_R2_allF = zeros(data.nChans,1);

data.studyRTCorr2mMagVec = zeros(data.nChans,data.nFreqs);
data.studyRTmVec_fastslow = zeros(data.nChans,data.nFreqs,2,2);
data.testRTCorr2mMagVec  = zeros(data.nChans,data.nFreqs);
data.testRTmVec_fastslow = zeros(data.nChans,data.nFreqs,2,2);

% correlation of quantities to RTs
trials = data.conds{13} ;
rts1   = -log10(data.studyRTs(trials));
rts2   = -log10(data.testRTs(trials));
trials_sFS{1} = trials &  (data.studyRTs<=data.studyRTsQuants(1));
trials_sFS{2} = trials &  (data.studyRTs>=data.studyRTsQuants(2));
trials_tFS{1} = trials & (data.testRTs<=data.testRTsQuants(1));
trials_tFS{2} = trials & (data.testRTs>=data.testRTsQuants(2));
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
        v1 = pca(cMI(trials,:)','numcomponents',5);
        m=fitglm(v1,rts1);       
        data.studyRTsMI_R2(cc,ff) =   m.Rsquared.Ordinary;
        m=fitglm(v1,rts2);       
        data.testRTsMI_R2(cc,ff)  =   m.Rsquared.Ordinary;
        
        % get mean trial direction 
        yy = mean(x.*exp(1j*phase),2);
        data.mMagVec(cc,ff,:) =  abs(yy);
        data.mPhVec(cc,ff,:)  =  angle(yy);
        
        % for trials of interest:
        data.studyRTCorr2mMagVec(cc,ff) = corr(rts1, abs(yy(trials)),'type','spearman');
        data.testRTCorr2mMagVec(cc,ff) = corr(rts2, abs(yy(trials)),'type','spearman');
        
        for ii = 1:2
            z=mean(yy(trials_sFS{ii}));
            data.studyRTmVec_fastslow(cc,ff,ii,:) = [abs(z) angle(z)];
            z=mean(yy(trials_tFS{ii}));
            data.testRTmVec_fastslow(cc,ff,ii,:)  = [abs(z) angle(z)];
        end
    end 
        X = squeeze(data.MI(cc,:,trials,:));
        X = permute(X,[2,1,3]); X = X(:,:);
        v1 = pca(X','numcomponents',12);        
        m=fitglm(v1,rts1);       
        data.studyRTsMI_R2_allF(cc) =   m.Rsquared.Ordinary;
        m=fitglm(v1,rts2);       
        data.testRTsMI_R2_allF(cc)  =   m.Rsquared.Ordinary;
        
        X= squeeze(data.mMagVec(cc,:,trials));
        v1 = pca(X,'numcomponents',12);        
        m=fitglm(v1,rts1);       
        data.studyRTsmMagVec_R2_allF(cc) =   m.Rsquared.Ordinary;
        m=fitglm(v1,rts2);       
        data.testRTsmMagVec_R2_allF(cc)  =   m.Rsquared.Ordinary;
end
return