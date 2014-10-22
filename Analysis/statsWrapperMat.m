function out = statsWrapperMat(X,groups)
%
% X -> matrix array
% the data matrix should have the same number of columns (time dimension)
% the rows represent channels
%
% the second level is on a channel grouping basis 
% groups is a vector indicating which channels are to be grouped, indicated by 
% a cardinal numnber
%
% group based statistics are baed on mean tests (i.e. t-tests)
%

nChans = size(X,1);
assert(nChans==numel(groups),'number of channels in groups does not match data')

nBins = size(X,2);

% get channel wise scores and p-values
out = [];
out.chanScore = X;

% get group statistcs
out.groups 			= groups;
groupIDs            = unique(groups);
out.nGroups 		= numel(groupIDs);
out.groupScores 	= zeros(out.nGroups,nBins);
out.groupPVals 		= zeros(out.nGroups,nBins);
for iGroup = 1:out.nGroups
	[~,p,~,t] = ttest(out.chanScore(groups==iGroup,:));
	out.groupScores(iGroup,:) 	= t.tstat;
	out.groupPVals(iGroup,:) 	= p;
end

% get two-way interactions 
out.GroupCombs 		= combnk(1:out.nGroups,2); 
out.nCombs 			= size(out.GroupCombs,1);
out.combScores 		= zeros(out.nCombs,nBins);
out.combPVals		= zeros(out.nCombs,nBins);
for iCombs = 1:out.nCombs
	comb = out.GroupCombs(iCombs,:);
	[~,p,~,t] = ttest2(out.chanScore(groups==comb(1),:),out.chanScore(groups==comb(2),:));
	out.combScores(iCombs,:) 	= t.tstat;
	out.combPVals(iCombs,:) 	= p;
end

return