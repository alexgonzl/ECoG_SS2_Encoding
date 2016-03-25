function out = PCAdecomp(data,opts)
% function that decomposes spectro temporal channel data.
% it returns the components, mapping of the componenets to channels
% variance explained
%
% 	opts:
% 		nBoot 		-> boot
%		KFold 	  	-> number of XVal folds.
% 		RTset 		-> {'studyRTs','testRTs'}