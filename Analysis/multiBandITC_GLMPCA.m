function out = multiBandITC_GLMPCA(data,opts)



out = [];
out.NumComponents = opts.NumComponents;
copyFields = {'nChans','SR','nEvents','trialTime','trialDur', ...
			'nTrialSamps','nFreqs','behavior','rois','conds','studyRTs',...
			'testRTs'};
for ii = 1:numel(copyFields)
	out.(copyFields{ii}) = data.(copyFields{ii});
end

trials = out.conds{13} ;
nTrials = sum(trials);
nChans = out.nChans;
nFreqs = out.nFreqs;

out.studyRTsPhaseGLMsPCA_R2     = zeros(nChans,nFreqs);
out.testRTsPhaseGLMsPCA_R2      = zeros(nChans,nFreqs);
out.PhaseGLMsPCAComps           = zeros(nChans,nFreqs,nTrials,opts.NumComponents);

out.studyRTsPhaseR2_allF        = zeros(nChans,1);
out.testRTsPhaseR2_allF         = zeros(nChans,1);
out.PhaseGLMsPCAComps_allF      = zeros(nChans,nTrials,opts.NumComponents);

rts1   = -log10(out.studyRTs(trials));
rts2   = -log10(out.testRTs(trials));

for ch = 1:nChans
    for ff = 1:nFreqs
        alpha = squeeze(data.trPh(ch,ff,trials,:));
        d_alpha = mod(diff(alpha,1,2),2*pi);
        v1=pca([cos(alpha) sin(alpha) d_alpha]','NumComponents',out.NumComponents);
        
        % study
        m=fitglm(v1,rts1);
        out.studyRTsPhaseGLMsPCA_R2(ch,ff) = m.Rsquared.Ordinary;

        % test
        m=fitglm(v1,rts2);
        out.testRTsPhaseGLMsPCA_R2(ch,ff) = m.Rsquared.Ordinary;
        out.studyRTsPhaseGLMsPCA_Comps(ch,ff,:,:) = v1;
        out.PhaseGLMsPCA_Comps(ch,ff,:,:) = v1;
    end
    alpha = squeeze(data.trPh(ch,:,trials,:));
    alpha2=permute(alpha,[2 1 3]); alpha2=alpha2(:,:);
    d_alpha = mod(diff(alpha2,1,2),2*pi);
    v1=pca([cos(alpha2) sin(alpha2) d_alpha]','NumComponents',out.NumComponents);
    m=fitglm(v1,rts1);
    out.studyRTsPhaseR2_allF(ch) = m.Rsquared.Ordinary;
    m=fitglm(v1,rts2);
    out.testRTsPhaseR2_allF(ch) = m.Rsquared.Ordinary;
    
    out.PhaseGLMsPCAComps_allF(ch,:,:) = v1;
    %fprintf(' Analyses for chan %i completed \n',ch)
end