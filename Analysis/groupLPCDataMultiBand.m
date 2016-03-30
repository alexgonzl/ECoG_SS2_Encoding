function data = groupLPCDataMultiBand(opts)

nSubjs = 7;
bands  = opts.bands;

extension = [opts.lockType 'sublogPower' opts.reference];

% store basic info.
fileName    = [opts.hems 'ERSPs' bands{1} 'Group' extension];
bandDat     = load([opts.dataPath fileName]);
data = [];
data.extension = extension;
fields = {'subjChans','LPCchanId','hemChan','ROIid','ROIs','nSubjLPCchans',...
    'behavior','testRTs','studyRTs','testRTs','Bins','trialDur','winSize','nBins','nChans',...
    'nTrialsInAnalysis','dataToStudyRTsCorr','dataToStudyRTsCorrP','dataToTestRTsCorr','dataToTestRTsCorrP'};
for ff = 1:numel(fields)
    data.(fields{ff})= bandDat.data.(fields{ff});
end
data.conds = bandDat.data.conds;
data.nBands = numel(bands);

data.BinERPs                = cell(nSubjs,1);
data.mBinChResp             = zeros(data.nBands,data.nChans,data.nBins); 
data.tBinChResp             = zeros(data.nBands,data.nChans,data.nBins); 
data.chBinScore             = zeros(data.nBands,data.nChans);
data.dataToStudyRTsCorr     = zeros(data.nBands,data.nChans,data.nBins);
data.dataToStudyRTsCorrP    = zeros(data.nBands,data.nChans,data.nBins);
data.dataToTestRTsCorr      = zeros(data.nBands,data.nChans,data.nBins);
data.dataToTestRTsCorrP     = zeros(data.nBands,data.nChans,data.nBins);
data.AnalysisBins           = bandDat.data.AnalysisBins;

for ba = 1:data.nBands
    fileName    = [opts.hems 'ERSPs' bands{ba} 'Group' extension];
    bandDat=load([opts.dataPath fileName]);
    for ss = 1:nSubjs
        trials = bandDat.data.conds{ss,13};
        data.trials{ss} = trials;
        data.BinERPs{ss}(ba,:,:,:) = bandDat.data.BinERP{ss}(:,trials,:);                
    end
    data.mBinChResp(ba,:,:)             = bandDat.data.mBinChResp;           
    data.tBinChResp(ba,:,:)             = bandDat.data.tBinChResp;           
    data.chBinScore(ba,:)               = bandDat.data.chBinScore;           
    data.dataToStudyRTsCorr(ba,:,:)     = bandDat.data.dataToStudyRTsCorr;
    data.dataToStudyRTsCorrP(ba,:,:)    = bandDat.data.dataToStudyRTsCorrP;
    data.dataToTestRTsCorr(ba,:,:)      = bandDat.data.dataToTestRTsCorr;
    data.dataToTestRTsCorrP(ba,:,:)     = bandDat.data.dataToTestRTsCorrP;
end

data.mBinROI_tChResp  = zeros(3,data.nBands,data.nBins);
data.tBinROI_tChResp  = zeros(3,data.nBands,data.nBins);
for rr = 1:3
    chans = data.ROIid==rr;
    for ba= 1:data.nBands
        temp=squeeze(data.tBinChResp(ba,chans,:));
        data.mBinROI_tChResp(rr,ba,:)=mean(temp);
        [~,~,~,t] = ttest(temp);
        data.tBinROI_tChResp(rr,ba,:)=t.tstat;
    end
end

% PCA analysis by trial. 
data.PCAtrialDecomp                 = [];
data.PCAtrialDecomp.nComps          = 10; nComps = data.PCAtrialDecomp.nComps;
data.PCAtrialDecomp.Comps           = cell(nSubjs,1);
data.PCAtrialDecomp.Projections     = cell(nSubjs,1);
data.PCAtrialDecomp.VarExp          = zeros(data.nChans,nComps);
data.PCAtrialDecomp.CorrStudyRTs    = zeros(data.nChans,nComps);
data.PCAtrialDecomp.CorrTestRTs     = zeros(data.nChans,nComps);

for ss=1:nSubjs
    x = permute(data.BinERPs{ss}(:,:,:,data.AnalysisBins),[2 3 1 4]);
    x = x(:,:,:);
    [d1,d2,d3] = size(x);
    rts1 = data.studyRTs{ss}(data.trials{ss});
    rts2 = data.testRTs{ss}(data.trials{ss});
    
    data.PCAtrialDecomp.Comps{ss} = zeros(d1,d2,nComps);
    data.PCAtrialDecomp.Projections{ss} = zeros(d1,d3,nComps);
    
    subjChans = find(data.subjChans==ss);
    for c=1:d1
        xx = squeeze(x(c,:,:));
        [C,S,E]= ppca(xx',nComps);  
        data.PCAtrialDecomp.Comps{ss}(c,:,:) = C;
        data.PCAtrialDecomp.Projections{ss}(c,:,:) = S;
        data.PCAtrialDecomp.VarExp(subjChans(c),:) = E;
        
        data.PCAtrialDecomp.CorrStudyRTs(subjChans(c),:) = corr(C,rts1);
        data.PCAtrialDecomp.CorrTestRTs(subjChans(c),:)  = corr(C,rts2);

    end
end

