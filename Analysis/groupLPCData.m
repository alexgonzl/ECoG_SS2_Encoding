function data = groupLPCData(opts)
% groups LPC channel data by subject.
% takes outputs from calcERP or calcERSP
%
% dependencies:
%       getBinSamps
%       signtest2
%       ranksum2

% data parameters
subjects    = opts.subjects;
reference   = opts.reference;
lockType    = opts.lockType;
band        = opts.band;

% erps or a spectral band
switch band
    case 'erp'
        baselineType = 'sub';
        analysisType = 'Amp';
        type    = 'ERP_Data/';
        preFix = 'ERPs' ;
    otherwise
        baselineType = 'sub';
        analysisType = 'logPower';
        type    = ['Spectral_Data/' band '/'];
        preFix = ['ERSPs' band];
end

extension = [lockType 'Lock' baselineType analysisType reference] ;
fileName = [preFix extension];

nSubjs  = numel(subjects);

data                = [];
data.options        = opts;
data.prefix         = preFix;
data.extension      = extension;
data.options        = opts;
data.subjChans      = [];
data.ERP            = [];
data.LPCchanId      = [];
data.ROIid          = []; data.ROIs = {'IPS','SPL','AG'};
data.subROIid       = []; data.subROIs = {'pIPS','aIPS','pSPL','aSPL'};
data.RTquantiles    = [0.3 0.7]; % quanntiles for analysis

for s = 1:nSubjs
    
    dataIn = load([opts.dataPath subjects{s} '/' type fileName '.mat']);
    
    temp = subjChanInfo(subjects{s});
    data.nSubjLPCchans(s)   = numel(temp.LPC);
    data.LPCchanId          = [data.LPCchanId; temp.LPC'];
    data.subjChans          = [data.subjChans s*ones(1,data.nSubjLPCchans(s))];
    
    ROIid = 1*ismember(temp.LPC,temp.IPS) + ...
        2*ismember(temp.LPC,temp.SPL) + ...
        3*ismember(temp.LPC,temp.AG);
    subROIid = 1*ismember(temp.LPC,temp.pIPS) + ...
        2*ismember(temp.LPC,temp.aIPS) + ...
        3*ismember(temp.LPC,temp.pSPL) + ...
        4*ismember(temp.LPC,temp.aSPL);
    
    data.ROIid              = [data.ROIid; ROIid'];
    data.subROIid           = [data.subROIid; subROIid'];
    
    data.ERP{s}             = squeeze(dataIn.data.erp(temp.LPC,:,:));
    
    data.remTrials{s}       = dataIn.data.behavior.subRem==1; % correctly remember trials
    data.testRTs{s}         = dataIn.data.behavior.testRTs;
    data.studyRTs{s}        = dataIn.data.behavior.studyRTs;
    data.testRTsQuants{s}   = quantile(data.testRTs{s}(data.remTrials{s}),data.RTquantiles);
    
    data.subFastRem{s}      = data.remTrials{s} & (data.testRTs{s} <= data.testRTsQuants{s}(1));
    data.subSlowRem{s}      = data.remTrials{s} & (data.testRTs{s} >= data.testRTsQuants{s}(2));
    
end

nChans = numel(data.subjChans);

% trial time info
data.baselineType   = dataIn.data.baselineType;
data.SR             = dataIn.data.SR;
data.nTrialSamps    =numel(dataIn.data.trialTime);
data.trialDur       = dataIn.data.trialDur;
data.trialTime      = linspace(data.trialDur(1),data.trialDur(2),data.nTrialSamps);

% bin info
data.winSize        = 0.1; % in seconds
data.sldWin         = data.winSize;
[data.BinSamps , data.Bins] = getBinSamps(data.winSize,data.sldWin,data.trialTime);
data.nBins          = size(data.Bins,1);

% pre allocation
fields = {'ZStat','PValZ','mCond1','mCond2','zCond1','zCond2','cCond1','cCond2',...
    'ZcCond1','ZcCond2','ZcStat'};

for f = fields
    data.(f{1}) = zeros(nChans,data.nTrialSamps);
end

fields2 = strcat('Bin',fields);
for f = fields2
    data.(f{1}) = zeros(nChans,data.nBins);
end

ch = 1;
data.BinERP     = cell(nSubjs,1);
for s = 1:nSubjs
    Cond1           = data.subFastRem{s};
    Cond2           = data.subSlowRem{s};
    RTs             = data.testRTs{s};
    
    nTrials         = numel(RTs);
    nSubjChans      = data.nSubjLPCchans(s);
    
    data.BinERP{s}      = nan(nSubjChans,nTrials,data.nBins);
    
    for Sch = 1: data.nSubjLPCchans(s)
        
        % original sampled data
        Z       = squeeze(data.ERP{s}(Sch,:,:));
        data    = getEffectScores(Z,RTs,Cond1,Cond2,data,ch,'');
        
        % binned
        binERP  = binTrials(Z,data.BinSamps);
        data.BinERP{s}(Sch,:,:) = binERP;
        data    = getEffectScores(binERP,RTs,Cond1,Cond2,data,ch,'Bin');
        
        ch = ch + 1;
    end
end

% channel by hemisphere id; 1 for lefts, 2 for rights
%data.hemChanId  = ismember(data.subjChans,find(strcmp(opts.hemId,'r')))'+1;

% % obtain group level roi main effects
% data.mainEfpValROIs = zeros(3,data.nBigBins,2);
% for hem = 1:2
%     for r   = 1: numel(data.ROIs)
%         chans = (data.ROIid == r) & (data.hemChanId == hem);
%         [~,data.mainEfpValROIs(r,:,hem)] = ttest(data.BinZStat(chans,:));
%     end
% end
%
% % obtain group level roi comparison statistics
% % three comparsions:
% data.ROIcontrasts       =  {'AG', 'IPS' ;'AG' ,'SPL';'IPS','SPL'};
% data.contrEfpValROIs    = zeros(3,data.nBins,2);
% for hem = 1:2
%     for rc = 1:size(data.ROIcontrasts,1)
%         r1 = find(strcmp(data.ROIcontrasts(rc,1),data.ROIs));
%         r2 = find(strcmp(data.ROIcontrasts(rc,2),data.ROIs));
%         chans1 = (data.ROIid == r1) & (data.hemChanId == hem);
%         chans2 = (data.ROIid == r2) & (data.hemChanId == hem);
%         [~,data.contrEfpValROIs(rc,:,hem)] = ttest2(data.BinZStat(chans1,:),data.BinZStat(chans2,:));
%     end
% end

% subject channel information
% for s = 1:nSubjs
%     subjChans = data.LPCchanId(data.subjChans==s);
%     temp = load(['./lib/elecLocs/subj' subjects{s} '_mni_elcoord_corrected.mat'],'mni_elcoord');
%     data.MNIlocsSubj{s} = temp.mni_elcoord(subjChans,:);
%     temp = load(['./lib/elecLocs/subj' subjects{s} '_electrodes_surface_loc_all1_correctnumbering.mat']);
%     data.OrigLocsSubj{s} = temp.elecmatrix(subjChans,:);
%     temp = load(['./lib/elecLocs/subj' subjects{s} '_cortex.mat']);
%     data.cortex{s} = temp.cortex;
% end
% 
% data.MNILocs = []; data.subjLocs = [];
% for s = 1:nSubjs
%     data.MNILocs = vertcat(data.MNILocs, data.MNIlocsSubj{s});
%     data.subjLocs = vertcat(data.subjLocs, data.OrigLocsSubj{s});
% end
% 
% switch data.options.hems
%     case 'l'
%         temp = load('./lib/elecLocs/MNI_cortex_left');
%         data.MNIcortex = temp.cortex;
%     case 'r'
%         temp = load('./lib/elecLocs/MNI_cortex_right');
%         data.MNIcortex = temp.cortex;
%     otherwise
%         temp = load('./lib/elecLocs/MNI_cortex_left');
%         data.lMNIcortex = temp.cortex;
%         temp = load('./lib/elecLocs/MNI_cortex_right');
%         data.rMNIcortex = temp.cortex;
%         temp = load('./lib/elecLocs/MNI_cortex_both');
%         data.MNIcortex = temp.cortex;
% end

end
function binZ = binTrials(Z,samps)
nT = size(Z,1);
nB = size(samps,1);

binZ = nan(nT,nB);
for tt = 1:nT
    for bb = 1:nB
        binZ(tt,bb) = mean(Z(tt,samps(bb,:)));
    end
end
end

function data = getEffectScores(Z,RTs,C1,C2,data,ch,preFix)
% data = getEffectScores(Z,RTs,H,CRs,data,preFix)
% inputs:
% Z     -> data matrix with 1st dimension of Ntrials. 2nd dim of timepoints
% RTs   -> vector with Ntrials
% C1    -> cond1, logical vector Ntrials
% C2   -> cond2, logical vector Ntrials
% data  -> data structure to write output to.
% preFix-> preFix string for the outputs

nC1    = sum(C1);
nC2    = sum(C2);

% original sampled data
X       = Z(C1,:);
Y       = Z(C2,:);

% store means for each condition
data.([preFix 'mCond1'])(ch,:)        = mean(X);
data.([preFix 'mCond2'])(ch,:)        = mean(Y);

if ~strcmp(preFix,'')
    % test each condition against 0
    [~,data.([preFix 'zCond1'])(ch,:)]                = signtest2(X);
    [~,data.([preFix 'zCond2'])(ch,:)]                = signtest2(Y);
    % test conditions against each other
    [data.([preFix 'PValZ'])(ch,:), data.([preFix 'ZStat'])(ch,:)] = ranksum2(X,Y);
    
    % correlate each condition to RT (transformed to z values)
    data.([preFix 'cCond1'])(ch,:)        = corr(X,log10(RTs(C1)));
    data.([preFix 'cCond2'])(ch,:)        = corr(Y,log10(RTs(C2)));
    data.([preFix 'ZcCond1'])(ch,:)       = atanh(data.([preFix 'cCond1'])(ch,:))*sqrt(nC1-3);
    data.([preFix 'ZcCond2'])(ch,:)       = atanh(data.([preFix 'cCond2'])(ch,:))*sqrt(nC2-3);
    data.([preFix 'ZcStat'])(ch,:)        = data.([preFix 'ZcCond1'])(ch,:)-data.([preFix 'ZcCond2'])(ch,:);
end
end