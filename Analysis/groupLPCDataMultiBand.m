function data = groupLPCDataMultiBand(opts)

nSubjs = 8;
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

data.avgDataTB_toStudyRTsCorr     = zeros(data.nBands,data.nChans,7);
data.avgDataTB_toStudyRTsCorrP    = zeros(data.nBands,data.nChans,7);
data.avgDataTB_toTestRTsCorr      = zeros(data.nBands,data.nChans,7);
data.avgDataTB_toTestRTsCorrP     = zeros(data.nBands,data.nChans,7);

data.AnalysisBins           = bandDat.data.AnalysisBins;
data.AnalysisBins2           = bandDat.data.AnalysisBins2;
data.Bins2                  = bandDat.data.Bins2;
for ba = 1:data.nBands
    fileName    = [opts.hems 'ERSPs' bands{ba} 'Group' extension];
    bandDat=load([opts.dataPath fileName]);
    for ss = 1:nSubjs
        trials = bandDat.data.conds{ss,13};
        data.trials{ss} = trials;
        data.BinERPs{ss}(ba,:,:,:) = bandDat.data.BinERP{ss}(:,trials,:);
        data.BinERPs2{ss}(ba,:,:,:) = bandDat.data.BinERP2{ss}(:,trials,:);
    end
    data.mBinChResp(ba,:,:)             = bandDat.data.mBinChResp;
    data.tBinChResp(ba,:,:)             = bandDat.data.tBinChResp;
    data.chBinScore(ba,:)               = bandDat.data.chBinScore;
    
    data.dataToStudyRTsCorr(ba,:,:)     = bandDat.data.dataToStudyRTsCorr;
    data.dataToStudyRTsCorrP(ba,:,:)    = bandDat.data.dataToStudyRTsCorrP;
    data.dataToTestRTsCorr(ba,:,:)      = bandDat.data.dataToTestRTsCorr;
    data.dataToTestRTsCorrP(ba,:,:)     = bandDat.data.dataToTestRTsCorrP;
    
    data.avgDataTB_toStudyRTsCorr(ba,:,:)     = bandDat.data.avgDataTB_toStudyRTsCorr;
    data.avgDataTB_toStudyRTsCorrP(ba,:,:)    = bandDat.data.avgDataTB_toStudyRTsCorrP;
    data.avgDataTB_toTestRTsCorr(ba,:,:)      = bandDat.data.avgDataTB_toTestRTsCorr;
    data.avgDataTB_toTestRTsCorrP(ba,:,:)     = bandDat.data.avgDataTB_toTestRTsCorrP;
    
    data.avgDataTB_PreToPostCorr(ba,:)          = bandDat.data.avgDataTB_PreToPostCorr;
    
    % models
    % pre post post-pre
    data.StudyRTs_PrePostActModel_Ts1(ba,:,:) = bandDat.data.StudyRTs_PrePostActModel_Ts1;
    data.TestRTs_PrePostActModel_Ts1(ba,:,:) = bandDat.data.TestRTs_PrePostActModel_Ts1;
    data.StudyRTs_PrePostActModel_Ts2(ba,:,:) = bandDat.data.StudyRTs_PrePostActModel_Ts2;
    data.TestRTs_PrePostActModel_Ts2(ba,:,:) = bandDat.data.TestRTs_PrePostActModel_Ts2;
    data.StudyRTs_PrePostActModel_Ts3(ba,:,:) = bandDat.data.StudyRTs_PrePostActModel_Ts3;
    data.TestRTs_PrePostActModel_Ts3(ba,:,:) = bandDat.data.TestRTs_PrePostActModel_Ts3;
    data.StudyRTs_PrePostActModel_Ts4(ba,:,:) = bandDat.data.StudyRTs_PrePostActModel_Ts4;
    data.TestRTs_PrePostActModel_Ts4(ba,:,:) = bandDat.data.TestRTs_PrePostActModel_Ts4;
    
end

data.mBinROI_tChResp  = zeros(3,data.nBands,data.nBins);
data.tBinROI_tChResp  = zeros(3,data.nBands,data.nBins);
data.pBinROI_tChResp  = zeros(3,data.nBands,data.nBins);
for rr = 1:3
    chans = data.ROIid==rr;
    for ba= 1:data.nBands
        temp=squeeze(data.tBinChResp(ba,chans,:));
        data.mBinROI_tChResp(rr,ba,:)=mean(temp);
        [~,p,~,t] = ttest(temp);
        data.tBinROI_tChResp(rr,ba,:)=t.tstat;
        data.pBinROI_tChResp(rr,ba,:)=p;
    end
end

%% pre to post analyses
nSubjs = 8;
preBins = 3;
postBins = 4:6;

nP = 2; % predictions are to study/test RTs
nConds = 3;
ROIs = {'IPS','SPL','AG'};
nROIs=3;

T = nan(nROIs,nROIs,nConds,nP,nSubjs);
P = nan(nROIs,nROIs,nConds,nP,nSubjs);

Study = cell(nROIs,1);
Test = cell(nROIs,1);

R_Vec = [];
Hem_Vec = [];

ROI_Combs = [];
Subj_Vec  = [];
PrePostC  = [];
for ss = 1:nSubjs
    
    subjROIchans = data.ROIid(data.subjChans==ss);
    trials = data.trials{ss};
    
    RTs = [];
    RTs{1} = log10(data.studyRTs{ss}(trials));
    RTs{2} = log10(data.testRTs{ss}(trials));
    
    nTrials = sum(trials);
    nChans = data.nSubjLPCchans(ss);
    % conds: pre post pre to post (#3)
    % predicting: studyRTs/testRTs (#2)
    bivariateChan_Act_To_RTCorr = nan(nChans,nChans,nConds,nP,nSubjs);
    
    x=squeeze(mean(data.BinERPs2{ss}(:,:,:,preBins),4));
    y=squeeze(mean(data.BinERPs2{ss}(:,:,:,postBins),4));
    
    n2 = nChans^2;
    n2_2 = (n2-nChans)/2;
    mask = tril(true(nChans),-1);
    
    C = cell(3,1);
    C{1} = zeros(nTrials,n2_2);
    C{2} = zeros(nTrials,n2_2);
    C{3} = zeros(nTrials,n2);
    
    M = cell(nChans,nChans);
    for ii=1:nChans
        for jj=1:nChans
            M{ii,jj} = strcat(ROIs{subjROIchans(ii)},'-',ROIs{subjROIchans(jj)});
        end
    end
    for ii = 1:nTrials
        xx=x(:,:,ii);
        yy=y(:,:,ii);
        c1 = atanh(corr(xx));
        c2 = atanh(corr(yy));
        c3 = atanh(corr(xx,yy));
        C{1}(ii,:) = c1(mask)';
        C{2}(ii,:) = c2(mask)';
        C{3}(ii,:) = c3(:)';
    end
    
    for pp =1:nP % study and test RTs
        for cc = 1:nConds % number of conditions (pre/post pre to post)
            if cc==3
                temp=reshape(corr(C{cc},RTs{pp},'type','spearman'),[nChans,nChans]);
                if pp==1
                    R_Vec       = [R_Vec;repmat(subjROIchans,[nChans,1])];
                    ROI_Combs   = [ROI_Combs;M(:)];
                    Subj_Vec    = [Subj_Vec;ss*ones(nChans^2,1)];
                    Hem_Vec     = [Hem_Vec; repmat(data.hemChan(data.subjChans==ss),[nChans,1])];
                    x = mean(C{3},1);
                    PrePostC    = [PrePostC;x'];
                end
            else
                temp = nan(nChans);
                temp(mask) =corr(C{cc},RTs{pp},'type','spearman');
            end
            if pp==1
                Study{cc}   = [Study{cc};temp(:)];
            else
                Test{cc} = [Test{cc}; temp(:)];
            end
        end
    end
    
end
%%

PrePostC2 = cell(size(PrePostC));
PrePostC2((PrePostC>=0)) = {'pos'};
PrePostC2((PrePostC<0)) = {'neg'};

tbl = table(Study{1},Study{2},Study{3},Test{1},Test{2},Test{3},PrePostC,PrePostC2,...
    ROI_Combs, categorical(Hem_Vec), categorical(Subj_Vec),'VariableNames',...
    {'S1','S2','S3','T1','T2','T3','PrePostC','PrePostC2','ROI_Pairs','Hem','Subjs'});

tbl2 = tbl(abs(PrePostC)>0.05,:);


data.Chan2ChanPrePostAct2RT = [];
data.Chan2ChanPrePostAct2RT.table = tbl;
data.Chan2ChanPrePostAct2RT.table2 = tbl2;

fprintf('Model fit to study RTs:')
m=fitlme(tbl2,'S3~1+ROI_Pairs*PrePostC2+Hem+(1|Subjs)');
m.anova()
data.Chan2ChanPrePostAct2RT.studyModel = m;
data.Chan2ChanPrePostAct2RT.studyANOVA = m.anova();

fprintf('Model fit to test RTs:')
m=fitlme(tbl2,'T3~1+ROI_Pairs*PrePostC2+Hem+(1|Subjs)');
m.anova()
data.Chan2ChanPrePostAct2RT.testModel = m;
data.Chan2ChanPrePostAct2RT.testANOVA = m.anova();

fprintf('Model fit to PrePostC:')
m=fitlme(tbl,'PrePostC~1+ROI_Pairs+Hem+(1|Subjs)');
m.anova()

data.Chan2ChanPrePostAct2RT.prepostCModel = m;
data.Chan2ChanPrePostAct2RT.prepostC_ANOVA = m.anova();

end
