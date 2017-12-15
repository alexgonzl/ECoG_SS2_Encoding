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

extension = [lockType baselineType analysisType reference] ;
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
data.hemChan        = [];
data.ROIid          = []; data.ROIs = {'IPS','SPL','AG','TPJ','SMG'};
data.subROIid       = []; data.subROIs = {'pIPS','aIPS','pSPL','aSPL'};
data.RTquantiles    = 0.1:0.1:0.9; % quanntiles for analysis

data.condNames      = {'abstract','concrete','absResp','conResp',...
    'oldResp','newResp', 'Rem','Forg', 'correctSemantic', 'correctStudyTest',...
    'studyRTs','testRTs','validEncRemem'};

for s = 1:nSubjs

    % load subjects data
    dataIn = load([opts.dataPath subjects{s} '/' type fileName '.mat']);

    % select, ID and store relevant ROI channels
    temp                    = subjChanInfo(subjects{s});
    data.subjHemsIds{s}        = temp.hemisphere;
    data.nSubjLPCchans(s)   = numel(temp.LPC);
    data.LPCchanId          = [data.LPCchanId; temp.LPC'];
    data.subjChans          = [data.subjChans; s*ones(data.nSubjLPCchans(s),1)];
    data.hemChan            = [data.hemChan; (strcmp(data.subjHemsIds{s},'r')+1)*ones(data.nSubjLPCchans(s),1)]; % one for lefts

    ROIid = 1*ismember(temp.LPC,temp.IPS) + ...
        2*ismember(temp.LPC,temp.SPL) + ...
        3*ismember(temp.LPC,temp.AG) + ...
        4*ismember(temp.LPC,temp.SMG) + ...
        5*ismember(temp.LPC,temp.TPJ) ;

    subROIid = 1*ismember(temp.LPC,temp.pIPS) + ...
        2*ismember(temp.LPC,temp.aIPS) + ...
        3*ismember(temp.LPC,temp.pSPL) + ...
        4*ismember(temp.LPC,temp.aSPL);

    data.ROIid              = [data.ROIid; ROIid'];
    data.subROIid           = [data.subROIid; subROIid'];

    % get data for LPC channels
    data.ERP{s}             = squeeze(dataIn.data.erp(temp.LPC,:,:));

    % behavior and conditions
    data.behavior{s}        = dataIn.data.behavior;
    data.conds{s,1}         = data.behavior{s}.cond==1  ;  % objective abstract
    data.conds{s,2}         = data.behavior{s}.cond==2  ;  % objective concrete
    data.conds{s,3}         = data.behavior{s}.sResp==1 ;  % abstract response
    data.conds{s,4}         = data.behavior{s}.sResp==2 ;  % concrete response
    data.conds{s,5}         = data.behavior{s}.tResp<=2 ;  % old response at test
    data.conds{s,6}         = data.behavior{s}.tResp>=2 ;  % new response at test
    data.conds{s,7}         = data.behavior{s}.subRem==1; % remembered
    data.conds{s,8}         = data.behavior{s}.subForg==1; % forgotten
    data.conds{s,9}         = (data.conds{s,1} & data.conds{s,3}) | ...
        (data.conds{s,2} & data.conds{s,4});              % correct semantic desc
    data.conds{s,10}        = data.conds{s,9} & data.conds{s,7}; % correct at study and test
    data.conds{s,11}        = data.behavior{s}.studyRTs ;
    data.conds{s,12}        = data.behavior{s}.testRTs ;
    data.conds{s,13}        = data.conds{s,7}&~isnan(data.conds{s,11}); % valid encoding, remembered trials

    data.studyIDatTest{s}   = data.behavior{s}.studyIDatTest;
    data.testRTs{s}         = dataIn.data.behavior.testRTs;
    data.studyRTs{s}        = dataIn.data.behavior.studyRTs;
    data.testRTsQuants{s}   = quantile(data.testRTs{s}(data.conds{s,7}==1),data.RTquantiles);

end

% trial time info
data.baselineType   = dataIn.data.baselineType;
data.SR             = dataIn.data.SR;
data.nTrialSamps    =numel(dataIn.data.trialTime);
data.trialDur       = dataIn.data.trialDur;
data.trialTime      = linspace(data.trialDur(1),data.trialDur(2),data.nTrialSamps);
% bin info
data.winSize                = 0.05; % in seconds
data.sldWin                 = data.winSize;
data.winSize2                = 0.50; % in seconds
data.sldWin2                 = data.winSize2;

% all bins
[data.BinSamps , data.Bins] = getBinSamps(data.winSize,data.sldWin,data.trialTime);
data.nBins                  = size(data.Bins,1);
[data.BinSamps2 , data.Bins2] = getBinSamps(data.winSize2,data.sldWin2,data.trialTime);
data.nBins2                  = size(data.Bins2,1);

% Get samples to be used in analyses. Excludes baseline period.
% RT-locked data has baseline outside the epochs. preStim has no baseline.
% other cases should include the samples greater than the last sample
% included in the baseline
if any(strcmp(opts.lockType,{'RT','preStim'}))
    data.validAnalysisSamples = true(size(data.trialTime));
    % analysis bins
    data.AnalysisBins = true(size(data.Bins(:,2),1),1);
    data.AnalysisBins2 = true(size(data.Bins2(:,2),1),1);
else
    data.baseLine       = dataIn.data.baseLine;
    data.validAnalysisSamples = data.trialTime>data.baseLine(2);
    data.AnalysisBins         = data.Bins(:,2)>=data.baseLine(2);

end
data.nValidAnalysisSamps = sum(data.validAnalysisSamples);

data.PostBins = data.Bins(:,2)>0.01;
data.PostBins1 = data.Bins(:,2)>0.01 & data.Bins(:,1)<0.51;
data.PostBins2 = data.Bins(:,2)>0.51 & data.Bins(:,1)<1.01;
data.PreBins  = data.Bins(:,1)>-0.51 & data.Bins(:,2)<0.01;

% number of channels
data.nChans                 = sum(data.nSubjLPCchans);

% Pre-allocations
data.BinERP                 = cell(nSubjs,1);

% all these correspond to the condition of interest: #13: valid encoding
% trial that were remembered.
data.condOfInterest         = 13;
data.mChResp                = zeros(data.nChans,data.nTrialSamps);
data.tChResp                = zeros(data.nChans,data.nTrialSamps);
data.chScore                = zeros(data.nChans,1); % for valid samples
data.mBinChResp             = zeros(data.nChans,data.nBins);
data.tBinChResp             = zeros(data.nChans,data.nBins);
data.chBinScore             = zeros(data.nChans,1);

data.nTrialsInAnalysis      = zeros(nSubjs,1);
data.trialsInAnalysis       = cell(nSubjs,1);
data.dataToStudyRTsCorr     = zeros(data.nChans,data.nBins);
data.dataToStudyRTsCorrP    = zeros(data.nChans,data.nBins);
data.dataToTestRTsCorr      = zeros(data.nChans,data.nBins);
data.dataToTestRTsCorrP     = zeros(data.nChans,data.nBins);

data.avgDataTB_toRTsCorrs     = {'pre','post','post1','post2','post3','post-pre','all'};
nC = numel(data.avgDataTB_toRTsCorrs);
data.avgDataTB_toStudyRTsCorr    = zeros(data.nChans,nC);
data.avgDataTB_toStudyRTsCorrP    = zeros(data.nChans,nC);
data.avgDataTB_toTestRTsCorr    = zeros(data.nChans,nC);
data.avgDataTB_toTestRTsCorrP    = zeros(data.nChans,nC);

data.StudyRTs_PrePostActModel       = cell(data.nChans,4);
data.TestRTs_PrePostActModel       = cell(data.nChans,4);
data.StudyRTs_PrePostActModel_Ts1   = zeros(data.nChans,3);
data.TestRTs_PrePostActModel_Ts1    = zeros(data.nChans,3);
data.StudyRTs_PrePostActModel_Ts2   = zeros(data.nChans,3);
data.TestRTs_PrePostActModel_Ts2    = zeros(data.nChans,3);
data.StudyRTs_PrePostActModel_Ts3   = zeros(data.nChans,4);
data.TestRTs_PrePostActModel_Ts3    = zeros(data.nChans,4);
data.StudyRTs_PrePostActModel_Ts4   = zeros(data.nChans,2);
data.TestRTs_PrePostActModel_Ts4    = zeros(data.nChans,2);

ch = 1;
for s = 1:nSubjs
    % remembered items and correctly identified semantics
    trials = data.conds{s,data.condOfInterest };
    nChans = data.nSubjLPCchans(s);
    nTrials = numel(trials);
    data.nTrialsInAnalysis(s) = sum(trials);
    data.trialsInAnalysis{s} = trials;

    testRTs = data.testRTs{s}(trials);
    studyRTs = data.studyRTs{s}(trials);

    % pre-allocate
    data.BinERP{s}             = zeros(nChans,nTrials,data.nBins);
    data.BinERP2{s}             = zeros(nChans,nTrials,data.nBins2);

    for Sch = 1: data.nSubjLPCchans(s)
        Z       = squeeze(data.ERP{s}(Sch,:,:));
        ZB      = binTrials(Z,data.BinSamps);
        ZB2      = binTrials(Z,data.BinSamps2);
        data.BinERP{s}(Sch,:,:) = ZB;
        data.BinERP2{s}(Sch,:,:) = ZB2;
        Z = Z(trials,:); ZB = ZB(trials,:); ZB2 = ZB2(trials,:); 

        % summaries only for trials of interest
        % Continous data
        data.mChResp(ch,:)   = mean(Z);
        [~,~,~,t]=ttest(Z);
        data.tChResp(ch,:) = t.tstat;
        data.chScore(ch)   = norm(t.tstat(data.validAnalysisSamples),inf);

        % Bin data
        data.mBinChResp(ch,:) =mean(ZB);
        [~,~,~,t]=ttest(data.BinERP{s}(Sch,trials,:));
        data.tBinChResp(ch,:) = t.tstat;
        data.chBinScore(ch)   = norm(squeeze(t.tstat(data.AnalysisBins)),inf);

        % compute relationship of data to RTs
        [data.dataToStudyRTsCorr(ch,:),data.dataToStudyRTsCorrP(ch,:)] = ...
            corrDataToRTs(ZB,studyRTs);

        [data.dataToTestRTsCorr(ch,:),data.dataToTestRTsCorrP(ch,:)] = ...
            corrDataToRTs(ZB,testRTs);

        % data split by time bins
        all = mean(ZB2(:,[3:6]),2);
        post = mean(ZB2(:,[4 5 6]),2);
        pre = ZB2(:,3);
        post1 = ZB2(:,4);
        post2 = ZB2(:,5);
        post3 = ZB2(:,6);        
        post_pre = post-pre;
        post1_pre = post1-pre;

        tbl=array2table([pre post1 post2 post3 post post_pre post1_pre all log10(studyRTs) log10(testRTs)],...
            'VariableNames',{'pre','post1','post2', 'post3','post','post_pre',...
            'post1_pre','all', 'log_studyRTs','log_testRTs'});

        % correlations to pre/post activity
        [c,p]=corr([pre post post1 post2 post3 post_pre all], log10(studyRTs),'type','spearman');
        data.avgDataTB_toStudyRTsCorr(ch,:) = c;
        data.avgDataTB_toStudyRTsCorrP(ch,:) = p;

        [c,p]=corr([pre post post1 post2 post3 post_pre all], log10(testRTs),'type','spearman');
        data.avgDataTB_toTestRTsCorr(ch,:) = c;
        data.avgDataTB_toTestRTsCorrP(ch,:) = p;

        [c,p]=corr(post,pre,'type','spearman');

        data.avgDataTB_PreToPostCorr(ch) = c;
        data.avgDataTB_PreToPostCorrP(ch) = p;

        % modelfit for study RTs
        mdl=fitlm(tbl,'log_studyRTs~1+pre+post');
        data.StudyRTs_PrePostActModel{Sch,1} = mdl;
        data.StudyRTs_PrePostActModel_Ts1(ch,1)  = mdl.Coefficients.tStat('pre');
        data.StudyRTs_PrePostActModel_Ts1(ch,2)  = mdl.Coefficients.tStat('post');
        mdl=fitlm(tbl,'log_studyRTs~1+post_pre');
        data.StudyRTs_PrePostActModel_Ts1(ch,3)  = mdl.Coefficients.tStat('post_pre');

        % modelfit for test RTs
        mdl=fitlm(tbl,'log_testRTs~1+pre+post');
        data.TestRTs_PrePostActModel{Sch,1} = mdl;
        data.TestRTs_PrePostActModel_Ts1(ch,1)  = mdl.Coefficients.tStat('pre');
        data.TestRTs_PrePostActModel_Ts1(ch,2)  = mdl.Coefficients.tStat('post');
        mdl=fitlm(tbl,'log_testRTs~1+post_pre');
        data.TestRTs_PrePostActModel_Ts1(ch,3)  = mdl.Coefficients.tStat('post_pre');

        % modelfit for study RTs
        mdl=fitlm(tbl,'log_studyRTs~1+pre+post1');
        data.StudyRTs_PrePostActModel{Sch,2} = mdl;
        data.StudyRTs_PrePostActModel_Ts2(ch,1)  = mdl.Coefficients.tStat('pre');
        data.StudyRTs_PrePostActModel_Ts2(ch,2)  = mdl.Coefficients.tStat('post1');
        mdl=fitlm(tbl,'log_studyRTs~1+post1_pre');
        data.StudyRTs_PrePostActModel_Ts2(ch,3)  = mdl.Coefficients.tStat('post1_pre');

        % modelfit for test RTs
        mdl=fitlm(tbl,'log_testRTs~1+pre+post1');
        data.TestRTs_PrePostActModel{Sch,2} = mdl;
        data.TestRTs_PrePostActModel_Ts2(ch,1)  = mdl.Coefficients.tStat('pre');
        data.TestRTs_PrePostActModel_Ts2(ch,2)  = mdl.Coefficients.tStat('post1');
        mdl=fitlm(tbl,'log_testRTs~1+post1_pre');
        data.TestRTs_PrePostActModel_Ts2(ch,3)  = mdl.Coefficients.tStat('post1_pre');
        
        % model fit by time
        mdl=fitlm(tbl,'log_studyRTs~1+pre+post1+post2+post3');
        data.StudyRTs_PrePostActModel{Sch,3} = mdl;
        data.StudyRTs_PrePostActModel_Ts3(ch,1)  = mdl.Coefficients.tStat('pre');
        data.StudyRTs_PrePostActModel_Ts3(ch,2)  = mdl.Coefficients.tStat('post1');
        data.StudyRTs_PrePostActModel_Ts3(ch,3)  = mdl.Coefficients.tStat('post2');
        data.StudyRTs_PrePostActModel_Ts3(ch,4)  = mdl.Coefficients.tStat('post3');
                
        % modelfit for test RTs
        mdl=fitlm(tbl,'log_testRTs~1+pre+post1+post2+post3');
        data.TestRTs_PrePostActModel{Sch,3} = mdl;
        data.TestRTs_PrePostActModel_Ts3(ch,1)  = mdl.Coefficients.tStat('pre');
        data.TestRTs_PrePostActModel_Ts3(ch,2)  = mdl.Coefficients.tStat('post1');
        data.TestRTs_PrePostActModel_Ts3(ch,3)  = mdl.Coefficients.tStat('post2');       
        data.TestRTs_PrePostActModel_Ts3(ch,4)  = mdl.Coefficients.tStat('post3');       

        % model fit by all and post_pre
        mdl=fitlm(tbl,'log_studyRTs~1+all+post_pre');
        data.StudyRTs_PrePostActModel{Sch,4} = mdl;
        data.StudyRTs_PrePostActModel_Ts4(ch,1)  = mdl.Coefficients.tStat('all');
        data.StudyRTs_PrePostActModel_Ts4(ch,2)  = mdl.Coefficients.tStat('post_pre');
                        
        % modelfit for test RTs
        mdl=fitlm(tbl,'log_testRTs~1+all+post_pre');
        data.TestRTs_PrePostActModel{Sch,4} = mdl;
        data.TestRTs_PrePostActModel_Ts4(ch,1)  = mdl.Coefficients.tStat('all');
        data.TestRTs_PrePostActModel_Ts4(ch,2)  = mdl.Coefficients.tStat('post_pre');      

        % all parameters: pre post1 post2 post3 post all
        data.StudyRTs_PrePostActModel_Ts4(ch,1)  = mdl.Coefficients.tStat('all');
        data.StudyRTs_PrePostActModel_Ts4(ch,2)  = mdl.Coefficients.tStat('post_pre');
                        
        % modelfit for test RTs
        mdl=fitlm(tbl,'log_testRTs~1+all+post_pre');
        data.TestRTs_PrePostActModel{Sch,4} = mdl;
        data.TestRTs_PrePostActModel_Ts4(ch,1)  = mdl.Coefficients.tStat('all');
        data.TestRTs_PrePostActModel_Ts4(ch,2)  = mdl.Coefficients.tStat('post_pre'); 
        ch = ch + 1;
    end
end
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
