
function out = PCAdecomp(data,opts)
% function that decomposes spectro temporal channel data.
% it returns the components, mapping of the componenets to channels
% variance explained
%
% 	opts:
% 		nBoot 		-> boot
%		KFold 	  	-> number of XVal folds.
% 		analysis 	-> {'activity','studyRTs','testRTs'}

out = [];
out.opts  = opts;
out.chans = find(opts.chans);
out.bins  = data.AnalysisBins;
out.nBins = sum(out.bins);
switch opts.analysis
    case 'activity'
        X = permute(data.tBinChResp,[2 1 3]);        
    case 'studyRT'
        X=permute(data.dataToStudyRTsCorr,[2 1 3]);
    case 'testRT'
        X=permute(data.dataToTestRTsCorr,[2 1 3]);
end

X = X(out.chans,:,out.bins);
X = X(:,:);

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(X);