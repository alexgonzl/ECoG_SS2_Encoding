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




nBins = data.nBins;
mB = zeros(data.nChans,nBins);
FitCor_sRTs = zeros(data.nChans,nBins);
FitCor_tRTs = zeros(data.nChans,nBins);

B = cell(nSubjs,1);
D1 = cell(nSubjs,1);
D2 = cell(nSubjs,1);

nLags = 41;

Chan_XC= cell(nSubjs,1);
Chan_XC_sRTs= cell(nSubjs,1);
Chan_XC_tRTs= cell(nSubjs,1);
ROI_XC= cell(nSubjs,1);

ROI_XC_sRTs= nan(nSubjs,nROIs,nROIs,nLags);
ROI_XC_tRTs= nan(nSubjs,nROIs,nROIs,nLags);

tROI_XC_sRTs= nan(nSubjs,nROIs,nROIs,nLags);
tROI_XC_tRTs= nan(nSubjs,nROIs,nROIs,nLags);
 
chanCount = 1;
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
    
    %% quadratic fits to spectrograms.
    B{ss} = zeros(nChans,nTrials,nBins);    

    %fr = (1:6)'-3.5; fr2 = fr.^2; 
    fr = [[-1 -1 -1 1 1 1]',ones(6,1)];
    for ii = 1:nChans
        for tr = 1:nTrials         
            A = squeeze(data.BinERPs{ss}(:,ii,tr,:));
            temp =  fr\A;
            B{ss}(ii,tr,:) = temp(1,:);  
        end
        mB(chanCount,:) = mean(B{ss}(ii,:,:));
        FitCor_sRTs(chanCount,:) = corr(squeeze(B{ss}(ii,:,:)),RTs{1},'type','spearman');        
        FitCor_tRTs(chanCount,:) = corr(squeeze(B{ss}(ii,:,:)),RTs{2},'type','spearman');
        chanCount = chanCount+1;
        %B(ii,:,:) = linSpectrogramFit(squeeze(data.BinERPs{ss}(:,ii,:,:)));
        %B2(ii,:,:,:) = quadSpectrogramFit(squeeze(data.BinERPs{ss}(:,ii,:,:)));
    end
    %%
    D1{ss} = nan(nChans,nChans,nTrials);
    D2{ss} = nan(nChans,nChans,nTrials);
    Chan_XC{ss} = nan(nChans,nChans,nTrials,nLags);
    Chan_XC_sRTs{ss} = nan(nChans,nChans,nLags);
    Chan_XC_tRTs{ss} = nan(nChans,nChans,nLags);
    ROI_XC{ss} = nan(nROIs,nROIs,nTrials,nLags);  
    
    for ii = 1:nChans
        for jj = 1:nChans
            for tt = 1:nTrials
                 %[cc,lags] = crosscorr(squeeze(B(ii,tt,:)),squeeze(B(jj,tt,:)));
                 [cc,lags] = crosscorr(squeeze(B{ss}(ii,tt,:)),squeeze(B{ss}(jj,tt,:)));
                 Chan_XC{ss}(ii,jj,tt,:) = cc;
                 [D1{ss}(ii,jj,tt),mi] = max(cc);
                 D2{ss}(ii,jj,tt) = lags(mi);
            end
            % put in lag part here, to obtain max lag per chan pair.
            Chan_XC_sRTs{ss}(ii,jj,:) = corr(squeeze(Chan_XC{ss}(ii,jj,:,:)),RTs{1},'type','spearman');
            Chan_XC_tRTs{ss}(ii,jj,:) = corr(squeeze(Chan_XC{ss}(ii,jj,:,:)),RTs{2},'type','spearman');            
        end
        
    end
    %%
    for r1 = 1:nROIs
        for r2 = 1:nROIs
            x= Chan_XC{ss}(subjROIchans==r1,subjROIchans==r2,:,:);
            y = squeeze(mean(mean(x,1),2));
            ROI_XC{ss}(r1,r2,:,:) = y;            
            ROI_XC_sRTs(ss,r1,r2,:) = corr(y,RTs{1},'type','spearman');   
            ROI_XC_tRTs(ss,r1,r2,:) = corr(y,RTs{2},'type','spearman');            
            
            x = Chan_XC_sRTs{ss}(subjROIchans==r1,subjROIchans==r2,:);
            xx = permute(x,[3,1,2]); xx = xx(:,:)';
            [~,~,~,t] = ttest(xx);
            tROI_XC_sRTs(ss,r1,r2,:) = t.tstat;
            
            x = Chan_XC_tRTs{ss}(subjROIchans==r1,subjROIchans==r2,:);
            xx = permute(x,[3,1,2]); xx = xx(:,:)';
            [~,~,~,t] = ttest(xx);
            tROI_XC_tRTs(ss,r1,r2,:) = t.tstat;
        end
    end
    %x
    
end

data.SpecLinFits.Fits =B;
data.SpecLinFits.ChansMeans = mB;
data.SpecLinFits.FitCor_sRTs =FitCor_sRTs;
data.SpecLinFits.FitCor_tRTs = FitCor_tRTs;
data.SpecLinFits.Chan_XC_sRTs = Chan_XC_sRTs;
data.SpecLinFits.Chan_XC_tRTs = Chan_XC_tRTs;
data.SpecLinFits.ROI_XC = ROI_XC;
data.SpecLinFits.ROI_XC_sRTs = ROI_XC_sRTs;
data.SpecLinFits.ROI_XC_tRTs = ROI_XC_tRTs;

data.SpecLinFits.tROI_XC_sRTs = tROI_XC_sRTs;
data.SpecLinFits.tROI_XC_tRTs = tROI_XC_tRTs;

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

