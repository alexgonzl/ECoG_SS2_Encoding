function data=reReferenceData(data, reReference, nRefChans, CARch)
% function that reReferences the data.
% nRefChans gets overwritten to zero if the re referencing doesn't need
% that parameter.
%
% no file dependencies.
% data should be the output of preProcessRawData.m

data.reReferencing = reReference; %{'allChCAR','LPCChCAR','origCAR'}
%% rereference channels and save
refChan = data.chanInfo.refChannel; data.refChan = refChan;

% common averagening
switch data.reReferencing
    case 'LPCChCAR'
        CARch = data.chanInfo.LPC;
    case 'allChCAR'
        CARch = data.chanInfo.CARChannels;
    case 'nonLPCCh'        
        CARch = data.chanInfo.other;
    case 'nonLPCleastVarCh'
        % sort by variance
        data.chanVar = var(data.signal,0,2);
        [~,i]=sort(data.chanVar);
        temp = i(ismember(i,data.chanInfo.other));
        CARch = temp(1:nRefChans);
    case 'nLPClowEvokedVar'
        [~,i]=sort(median(data.evokedVar,2));
        temp = i(ismember(i,data.chanInfo.other));
        CARch = temp(1:nRefChans);
    case 'nonLPCleasL1TvalCh'
        [~,i]=sort(data.RefChanTotalTstat);
        temp = i(ismember(i,data.chanInfo.other))';
        CARch = temp(1:nRefChans);        
    case 'other' 
        % give specific channels
    otherwise 
        % don't do anything.
        CARch = refChan;
end

data.nRefChans = numel(CARch);
data.CARch = CARch;
data.CARSignal = mean(data.signal(CARch,:),1);
data.signal = bsxfun(@minus,data.signal,data.CARSignal);
data.signal(refChan,:) = data.CARSignal;

