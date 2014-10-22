function out = BiVariateStatsWrapper(X,Y,groups)
%
% X & Y -> cell array containing a data matrix per channel
% the data matrix should have the same number of columns (time dimension)
% the rows represent trials, which can be different by channel (different subjects)
%
% Analysis are first performed on a channel basis to obtain a channel wise summary 
%
% the second level is on a channel grouping basis 
% groups is a vector indicating which channels are to be grouped, indicated by 
% a cardinal numnber
%
% channel wise statistics are a based on median differences (i.e. ranksum/ wilcoxon tests)
% group based statistics are baed on mean tests (i.e. t-tests)
%

% basic assertions about the data:
assert(iscell(X) && iscell(X),'incorrect input type')

nChans = numel(X);
assert(nChans==numel(Y),'number of channels does not match between X and Y')
assert(nChans==numel(groups),'number of channels in groups does not match data')

nBins = size(X{1},2);
assert(nBins==size(Y{1},2),'number of channels does not match between X and Y')

% get channel wise scores and p-values
out = [];
out.chanScores = zeros(nChans,nBins);
out.chanPVals = zeros(nChans,nBins);

for iChans = 1:nChans	
    [out.chanPVals(iChans,:),out.chanScore(iChans,:)] = ranksum2(X{iChans},Y{iChans});
end

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