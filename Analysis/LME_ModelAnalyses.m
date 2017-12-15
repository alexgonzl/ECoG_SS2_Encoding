
function out = LME_ModelAnalyses(data)


%% models with time bins as independent preds.
if 1
    %% LME model analyses on pre/post activity to RT
    % Study
    Y1=data.StudyRTs_PrePostActModel_Ts1;
    
    % bands x channel x pre/post/inter
    [nBands,nChans,nConds] = size(Y1);
    N = numel(Y1);
    study = Y1(:);
    
    % Test
    test=data.TestRTs_PrePostActModel_Ts1(:);
    
    % set up of independent variables to match data:
    % frequency bands
    bands = repmat((1:nBands)',[nChans*nConds,1]);
    bands_IDs = {'delta','theta','alpha','beta','lgam','hgam'};
    bands_s = cell(N,1);
    for ii = 1:6
        bands_s(bands==ii)= bands_IDs(ii);
    end
    
    % conditions pre=1,post=2, inter=3
    conds = [ones(nBands*nChans,1);2*ones(nBands*nChans,1);3*ones(nBands*nChans,1)];
    
    conds_s = cell(N,1);
    conds_s(conds==1)={'pre'};
    conds_s(conds==2)={'post'};
    conds_s(conds==3)={'post-pre'};
    
    % channels
    temp=repmat(1:nChans,[nBands,1]);
    chans=repmat(temp(:),[nConds,1]);
    
    % channels by subject ID
    nSubjs=8;
    subjChans = zeros(N,1);
    for ii=1:nSubjs
        chanIDs = find(data.subjChans==ii);
        subjChans(ismember(chans,chanIDs))=ii;
    end
    
    % channels by hemisphere
    hemChans = zeros(N,1);
    for ii=1:2
        chanIDs = find(data.hemChan==ii);
        hemChans(ismember(chans,chanIDs))=ii;
    end
    
    % region: ips, spl, ag
    rois = zeros(N,1);
    rois_s = cell(N,1);
    rois_ids = {'IPS','SPL','AG'};
    for ii=1:3
        chanIDs = find(data.ROIid==ii);
        rois(ismember(chans,chanIDs))=ii;
        rois_s(ismember(chans,chanIDs))=rois_ids(ii);
    end
    
    tbl=table(study,test,bands_s, conds_s, ...
        rois_s, categorical(subjChans), categorical(hemChans), categorical(chans),...
        'variablenames',{'study','test','bands','conds','rois','subj','hem','chans'});
    
    % store!
    data.PrePostActRT = [];
    data.PrePostActRT.nBa_nCh_nCo = [nBands,nChans,nConds];
    data.PrePostActRT.vecBands = bands;
    data.PrePostActRT.vecConds = conds;
    data.PrePostActRT.vecConds_s = conds_s;
    data.PrePostActRT.vecROIs = rois;
    data.PrePostActRT.vecROIs_s = rois_s;
    data.PrePostActRT.vecHemChans = hemChans;
    data.PrePostActRT.vecSubjChans = subjChans;
    data.PrePostActRT.vecChanID   = chans;
    
    data.PrePostActRT.vecTestDat = test;
    data.PrePostActRT.vecStudyDat = study;
    data.PrePostActRT.datTable = tbl;
    
    data.PrePostActRT.studyModel = fitlme(tbl,'study~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
    data.PrePostActRT.studyANOVA = anova(data.PrePostActRT.studyModel);
    
    data.PrePostActRT.testModel = fitlme(tbl,'test~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
    data.PrePostActRT.testANOVA = anova(data.PrePostActRT.testModel);
    
    % correct FDR correction for MC across BOTH models
    x= [data.PrePostActRT.studyANOVA.pValue; data.PrePostActRT.testANOVA.pValue];
    [~,a] = fdr_bh(x);
    
    data.PrePostActRT.AcrossANOVAs_q = a;
    p = data.PrePostActRT.studyANOVA.pValue;
    data.PrePostActRT.studyANOVA_FDR = data.PrePostActRT.studyANOVA(p<=a,:);
    p = data.PrePostActRT.testANOVA.pValue;
    data.PrePostActRT.testANOVA_FDR = data.PrePostActRT.testANOVA(p<=a,:);
    
    %% models for post1 which is right after stim onset 500ms.
    study=data.StudyRTs_PrePostActModel_Ts2(:);
    test=data.TestRTs_PrePostActModel_Ts2(:);
    tbl=table(study,test,bands_s, conds_s, ...
        rois_s, categorical(subjChans), categorical(hemChans), categorical(chans),...
        'variablenames',{'study','test','bands','conds','rois','subj','hem','chans'});
    
    data.PrePostActRT.vecTestDat2 = test;
    data.PrePostActRT.vecStudyDat2 = study;
    data.PrePostActRT.datTable2 = tbl;
    
    data.PrePostActRT.studyModel2 = fitlme(tbl,'study~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
    data.PrePostActRT.studyANOVA2 = anova(data.PrePostActRT.studyModel2);
    
    data.PrePostActRT.testModel2 = fitlme(tbl,'test~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
    data.PrePostActRT.testANOVA2 = anova(data.PrePostActRT.testModel2);
    
    % correct FDR correction for MC across BOTH models
    x= [data.PrePostActRT.studyANOVA2.pValue; data.PrePostActRT.testANOVA2.pValue];
    [~,a] = fdr_bh(x);
    
    data.PrePostActRT.AcrossANOVAs_q2 = a;
    p = data.PrePostActRT.studyANOVA2.pValue;
    data.PrePostActRT.studyANOVA_FDR2 = data.PrePostActRT.studyANOVA2(p<=a,:);
    p = data.PrePostActRT.testANOVA2.pValue;
    data.PrePostActRT.testANOVA_FDR2 = data.PrePostActRT.testANOVA2(p<=a,:);
    
    %% model 3: pre +post1 +post2+post3
    
    study=data.StudyRTs_PrePostActModel_Ts3(:);
    test=data.TestRTs_PrePostActModel_Ts3(:);
    
    [nBands,nChans,nConds] = size(data.StudyRTs_PrePostActModel_Ts3);
    N = numel(study);
    % conditions pre=1,post=2, inter=3
    conds = [ones(nBands*nChans,1);2*ones(nBands*nChans,1);3*ones(nBands*nChans,1);4*ones(nBands*nChans,1)];
  
    conds_s = cell(N,1);
    conds_s(conds==1)={'pre'};
    conds_s(conds==2)={'post1'};
    conds_s(conds==3)={'post2'};
    conds_s(conds==4)={'post3'};

    % set up of independent variables to match data:
    % frequency bands
    bands = repmat((1:nBands)',[nChans*nConds,1]);
    bands_IDs = {'delta','theta','alpha','beta','lgam','hgam'};
    bands_s = cell(N,1);
    for ii = 1:6
        bands_s(bands==ii)= bands_IDs(ii);
    end
    
    % channels
    temp=repmat(1:nChans,[nBands,1]);
    chans=repmat(temp(:),[nConds,1]);
    
    % channels by subject ID
    nSubjs=8;
    subjChans = zeros(N,1);
    for ii=1:nSubjs
        chanIDs = find(data.subjChans==ii);
        subjChans(ismember(chans,chanIDs))=ii;
    end
    
    % channels by hemisphere
    hemChans = zeros(N,1);
    for ii=1:2
        chanIDs = find(data.hemChan==ii);
        hemChans(ismember(chans,chanIDs))=ii;
    end
    
    % region: ips, spl, ag
    rois = zeros(N,1);
    rois_s = cell(N,1);
    rois_ids = {'IPS','SPL','AG'};
    for ii=1:3
        chanIDs = find(data.ROIid==ii);
        rois(ismember(chans,chanIDs))=ii;
        rois_s(ismember(chans,chanIDs))=rois_ids(ii);
    end
    
    tbl=table(study,test,bands_s, conds, conds_s, ...
        rois_s, categorical(subjChans), categorical(hemChans), categorical(chans),...
        'variablenames',{'study','test','bands','conds_cont', 'conds','rois','subj','hem','chans'});
    
    data.PrePostActRT.vecTestDat3 = test;
    data.PrePostActRT.vecStudyDat3 = study;
    data.PrePostActRT.datTable3 = tbl;
    
    data.PrePostActRT.studyModel3 = fitlme(tbl,'study~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
    data.PrePostActRT.studyANOVA3 = anova(data.PrePostActRT.studyModel3);
    
    data.PrePostActRT.testModel3 = fitlme(tbl,'test~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
    data.PrePostActRT.testANOVA3 = anova(data.PrePostActRT.testModel3);
    
    % correct FDR correction for MC across BOTH models
    x= [data.PrePostActRT.studyANOVA3.pValue; data.PrePostActRT.testANOVA3.pValue];
    [~,a] = fdr_bh(x);
    
    data.PrePostActRT.AcrossANOVAs_q3 = a;
    p = data.PrePostActRT.studyANOVA3.pValue;
    data.PrePostActRT.studyANOVA_FDR3 = data.PrePostActRT.studyANOVA3(p<=a,:);
    p = data.PrePostActRT.testANOVA3.pValue;
    data.PrePostActRT.testANOVA_FDR3 = data.PrePostActRT.testANOVA3(p<=a,:);
    
    %% %% LME model analyses on pre/post only, without the difference.
    % Study
    Y1=data.StudyRTs_PrePostActModel_Ts1(:,:,1:2);
    
    % bands x channel x pre/post/inter
    [nBands,nChans,nConds] = size(Y1);
    N = numel(Y1);
    study = Y1(:);
    
    % Test
    test=data.TestRTs_PrePostActModel_Ts1(:,:,1:2);
    test=test(:);
    
    % set up of independent variables to match data:
    % frequency bands
    bands = repmat((1:nBands)',[nChans*nConds,1]);
    bands_IDs = {'delta','theta','alpha','beta','lgam','hgam'};
    bands_s = cell(N,1);
    for ii = 1:6
        bands_s(bands==ii)= bands_IDs(ii);
    end
    
    % conditions pre=1,post=2, inter=3
    conds = [ones(nBands*nChans,1);2*ones(nBands*nChans,1);3*ones(nBands*nChans,1)];
    
    conds_s = cell(N,1);
    conds_s(conds==1)={'pre'};
    conds_s(conds==2)={'post'};
    %conds_s(conds==3)={'post-pre'};
    
    % channels
    temp=repmat(1:nChans,[nBands,1]);
    chans=repmat(temp(:),[nConds,1]);
    
    % channels by subject ID
    nSubjs=8;
    subjChans = zeros(N,1);
    for ii=1:nSubjs
        chanIDs = find(data.subjChans==ii);
        subjChans(ismember(chans,chanIDs))=ii;
    end
    
    % channels by hemisphere
    hemChans = zeros(N,1);
    for ii=1:2
        chanIDs = find(data.hemChan==ii);
        hemChans(ismember(chans,chanIDs))=ii;
    end
    
    % region: ips, spl, ag
    rois = zeros(N,1);
    rois_s = cell(N,1);
    rois_ids = {'IPS','SPL','AG'};
    for ii=1:3
        chanIDs = find(data.ROIid==ii);
        rois(ismember(chans,chanIDs))=ii;
        rois_s(ismember(chans,chanIDs))=rois_ids(ii);
    end
    
    tbl=table(study,test,bands,bands_s, conds_s, ...
        rois_s, categorical(subjChans), categorical(hemChans), categorical(chans),...
        'variablenames',{'study','test','bands_c','bands','conds','rois','subj','hem','chans'});
    
    
    data.PrePostActRT.datTableX = tbl;
    
    data.PrePostActRT.studyModelX = fitlme(tbl,'study~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
    data.PrePostActRT.studyANOVAX = anova(data.PrePostActRT.studyModelX);
    
    data.PrePostActRT.testModelX = fitlme(tbl,'test~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
    data.PrePostActRT.testANOVAX = anova(data.PrePostActRT.testModelX);
    
    % correct FDR correction for MC across BOTH models
    x= [data.PrePostActRT.studyANOVAX.pValue; data.PrePostActRT.testANOVAX.pValue];
    [~,a] = fdr_bh(x);
    
    data.PrePostActRT.AcrossANOVAs_qX = a;
    p = data.PrePostActRT.studyANOVAX.pValue;
    data.PrePostActRT.studyANOVA_FDRX = data.PrePostActRT.studyANOVAX(p<=a,:);
    p = data.PrePostActRT.testANOVAX.pValue;
    data.PrePostActRT.testANOVA_FDRX = data.PrePostActRT.testANOVAX(p<=a,:);
    
    
    %% using time as independent predictors.
else
    % Study
    % get pre(500ms) and post(1500ms) regressors
    pre_S = data.StudyRTs_PrePostActModel_Ts1(:,:,1);
    pre_S = pre_S(:);
    post_S = data.StudyRTs_PrePostActModel_Ts1(:,:,2);
    post_S = post_S(:);
    post_pre_S = data.StudyRTs_PrePostActModel_Ts1(:,:,3);
    post_pre_S = post_pre_S(:);    
    
    % bands x channel x pre/post/inter
    nBands = data.nBands;
    nChans = data.nChans;
    N = numel(pre_S);
    
    % Test
    test=data.TestRTs_PrePostActModel_Ts1(:);
    pre_T = data.TestRTs_PrePostActModel_Ts1(:,:,1); pre_T = pre_T(:);
    post_T = data.TestRTs_PrePostActModel_Ts1(:,:,2); post_T = post_T(:);
    post_pre_T = data.TestRTs_PrePostActModel_Ts1(:,:,3); post_pre_T = post_pre_T(:);
    
    % set up of independent variables to match data:
    % frequency bands
    bands = repmat((1:nBands)',[nChans,1]);
    bands_IDs = {'delta','theta','alpha','beta','lgam','hgam'};
    bands_s = cell(N,1);
    for ii = 1:6
        bands_s(bands==ii)= bands_IDs(ii);
    end   
    
    % channels
    chans=repmat(1:nChans,[nBands,1]);
    chans = chans(:);
    
    % channels by subject ID
    nSubjs=8;
    subjChans = zeros(N,1);
    for ii=1:nSubjs
        chanIDs = find(data.subjChans==ii);
        subjChans(ismember(chans,chanIDs))=ii;
    end
    
    % channels by hemisphere
    hemChans = zeros(N,1);
    for ii=1:2
        chanIDs = find(data.hemChan==ii);
        hemChans(ismember(chans,chanIDs))=ii;
    end
    
    % region: ips, spl, ag
    rois = zeros(N,1);
    rois_s = cell(N,1);
    rois_ids = {'IPS','SPL','AG'};
    for ii=1:3
        chanIDs = find(data.ROIid==ii);
        rois(ismember(chans,chanIDs))=ii;
        rois_s(ismember(chans,chanIDs))=rois_ids(ii);
    end
    
    tbl=table(pre_S,pre_T,post_S,post_T,post_pre_S,post_pre_T, bands, bands_s,...
        rois_s, categorical(subjChans), categorical(hemChans), categorical(chans),...
        'variablenames',{'pre_S','pre_T','post_S','post_T','post_pre_S','post_pre_T',...
        'bands_c','bands','rois','subj','hem','chans'});
    
    
    % store!
    data.PrePostActRT = [];
    data.PrePostActRT.nBa_nCh_nCo = [nBands,nChans,nConds];
    data.PrePostActRT.vecBands = bands;
    data.PrePostActRT.vecConds = conds;
    data.PrePostActRT.vecConds_s = conds_s;
    data.PrePostActRT.vecROIs = rois;
    data.PrePostActRT.vecROIs_s = rois_s;
    data.PrePostActRT.vecHemChans = hemChans;
    data.PrePostActRT.vecSubjChans = subjChans;
    data.PrePostActRT.vecChanID   = chans;
    
    data.PrePostActRT.vecTestDat = test;
    data.PrePostActRT.vecStudyDat = study;
    data.PrePostActRT.datTable = tbl;
    
    data.PrePostActRT.studyModel = fitlme(tbl,'study~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
    data.PrePostActRT.studyANOVA = anova(data.PrePostActRT.studyModel);
    
    data.PrePostActRT.testModel = fitlme(tbl,'test~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
    data.PrePostActRT.testANOVA = anova(data.PrePostActRT.testModel);
    
    % correct FDR correction for MC across BOTH models
    x= [data.PrePostActRT.studyANOVA.pValue; data.PrePostActRT.testANOVA.pValue];
    [~,a] = fdr_bh(x);
    
    data.PrePostActRT.AcrossANOVAs_q = a;
    p = data.PrePostActRT.studyANOVA.pValue;
    data.PrePostActRT.studyANOVA_FDR = data.PrePostActRT.studyANOVA(p<=a,:);
    p = data.PrePostActRT.testANOVA.pValue;
    data.PrePostActRT.testANOVA_FDR = data.PrePostActRT.testANOVA(p<=a,:);
end

%% 


out=data;

% %% models for post1 which is right after stim onset 500ms.
% study=data.StudyRTs_PrePostActModel_Ts2(:);
% test=data.TestRTs_PrePostActModel_Ts2(:);
% tbl=table(study,test,bands_s, conds_s, ...
%     rois_s, categorical(subjChans), categorical(hemChans), categorical(chans),...
%     'variablenames',{'study','test','bands','conds','rois','subj','hem','chans'});
% 
% data.PrePostActRT.vecTestDat2 = test;
% data.PrePostActRT.vecStudyDat2 = study;
% data.PrePostActRT.datTable2 = tbl;
% 
% data.PrePostActRT.studyModel2 = fitlme(tbl,'study~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
% data.PrePostActRT.studyANOVA2 = anova(data.PrePostActRT.studyModel2);
% 
% data.PrePostActRT.testModel2 = fitlme(tbl,'test~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
% data.PrePostActRT.testANOVA2 = anova(data.PrePostActRT.testModel2);
% 
% % correct FDR correction for MC across BOTH models
% x= [data.PrePostActRT.studyANOVA2.pValue; data.PrePostActRT.testANOVA2.pValue];
% [~,a] = fdr_bh(x);
% 
% data.PrePostActRT.AcrossANOVAs_q2 = a;
% p = data.PrePostActRT.studyANOVA2.pValue;
% data.PrePostActRT.studyANOVA_FDR2 = data.PrePostActRT.studyANOVA2(p<=a,:);
% p = data.PrePostActRT.testANOVA2.pValue;
% data.PrePostActRT.testANOVA_FDR2 = data.PrePostActRT.testANOVA2(p<=a,:);
% 
% %% model 3: pre +post1 +post2+post3
% 
% study=data.StudyRTs_PrePostActModel_Ts3(:);
% test=data.TestRTs_PrePostActModel_Ts3(:);
% conds_s = cell(N,1);
% conds_s(conds==1)={'pre'};
% conds_s(conds==2)={'post1'};
% conds_s(conds==3)={'post2'};
% 
% tbl=table(study,test,bands_s, conds_s, ...
%     rois_s, categorical(subjChans), categorical(hemChans), categorical(chans),...
%     'variablenames',{'study','test','bands','conds','rois','subj','hem','chans'});
% 
% data.PrePostActRT.vecTestDat3 = test;
% data.PrePostActRT.vecStudyDat3 = study;
% data.PrePostActRT.datTable3 = tbl;
% 
% data.PrePostActRT.studyModel3 = fitlme(tbl,'study~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
% data.PrePostActRT.studyANOVA3 = anova(data.PrePostActRT.studyModel3);
% 
% data.PrePostActRT.testModel3 = fitlme(tbl,'test~1+bands*conds*rois+hem+(1|chans)+(1|subj)');
% data.PrePostActRT.testANOVA3 = anova(data.PrePostActRT.testModel3);
% 
% % correct FDR correction for MC across BOTH models
% x= [data.PrePostActRT.studyANOVA3.pValue; data.PrePostActRT.testANOVA3.pValue];
% [~,a] = fdr_bh(x);
% 
% data.PrePostActRT.AcrossANOVAs_q3 = a;
% p = data.PrePostActRT.studyANOVA3.pValue;
% data.PrePostActRT.studyANOVA_FDR3 = data.PrePostActRT.studyANOVA3(p<=a,:);
% p = data.PrePostActRT.testANOVA3.pValue;
% data.PrePostActRT.testANOVA_FDR3 = data.PrePostActRT.testANOVA3(p<=a,:);


end

