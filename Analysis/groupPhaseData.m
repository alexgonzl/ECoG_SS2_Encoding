function data = groupPhaseData(opts)

subjects    = opts.subjects;
reference   = opts.reference;
nFreqs      = opts.nFreq;
lockType    = opts.lockType;
nSubjs      = numel(subjects);

fileName1 = ['multiBandITC_N' num2str(nFreqs) '_' lockType reference '.mat'];
fileName2 = ['multiBandITC_GLMPCA_N' num2str(nFreqs) '_' lockType reference '.mat'];
fileName3 = ['modIndex_' lockType reference '.mat'];

data                = [];
data.options        = opts;
data.subjChans      = [];
data.LPCchanId      = [];
data.hemChan        = [];
data.ROIid          = []; data.ROIs = {'IPS','SPL','AG','TPJ','SMG'};

data.condNames      = {'abstract','concrete','absResp','conResp',...
    'oldResp','newResp', 'Rem','Forg', 'correctSemantic', 'correctStudyTest',...
    'studyRTs','testRTs','validEncRemem'};

% allocate fields
itcFields = {'ITC','ITC_studyRTsplit','ITC_testRTsplit','studyRTsToPhaseCorr',...
    'testRTsToPhaseCorr'};

itcGLMFields = {'studyRTsPhaseGLMsPCA_R2','testRTsPhaseGLMsPCA_R2',...
    'studyRTsPhaseR2_allF','testRTsPhaseR2_allF'};

miFields    = {'studyRTsMI_R2','studyRTsMI_R2_allF','testRTsMI_R2','testRTsMI_R2_allF', ...
            'studyRTsmMagVec_R2_allF','testRTsmMagVec_R2_allF', ...
        'studyRTCorr2mMagVec','studyRTmVec_fastslow','testRTCorr2mMagVec','testRTmVec_fastslow'};
for ff = 1:numel(itcFields)
    data.(itcFields{ff}) = [];
end

for ff = 1:numel(itcGLMFields)
    data.(itcGLMFields{ff}) = [];
end

for ff = 1:numel(miFields)
    data.(miFields{ff}) = [];
end

for ss = 1:nSubjs
    % load subjects data
    dataIn = load([opts.dataPath subjects{ss} '/ITC_Data/' fileName1 ]);
    
    % select, ID and store relevant ROI channels
    temp                    = subjChanInfo(subjects{ss});
    data.subjHemsIds{ss}        = temp.hemisphere;
    data.nSubjLPCchans(ss)   = numel(temp.LPC);
    data.LPCchanId          = [data.LPCchanId; temp.LPC'];
    data.subjChans          = [data.subjChans; ss*ones(data.nSubjLPCchans(ss),1)];
    data.hemChan            = [data.hemChan; (strcmp(data.subjHemsIds{ss},'r')+1)*ones(data.nSubjLPCchans(ss),1)]; % one for lefts
    
    ROIid = 1*ismember(temp.LPC,temp.IPS) + ...
        2*ismember(temp.LPC,temp.SPL) + ...
        3*ismember(temp.LPC,temp.AG) + ...
        4*ismember(temp.LPC,temp.SMG) + ...
        5*ismember(temp.LPC,temp.TPJ) ;
    
    data.ROIid          = [data.ROIid; ROIid'];
    
    data.behavior       = dataIn.data.behavior;
    for co = 1:numel(data.condNames)
        data.conds{ss,co} = dataIn.data.conds{co};
    end
    data.testRTs{ss}          = dataIn.data.testRTs;
    data.studyRTs{ss}         = dataIn.data.studyRTs;
    data.studyRTsQuants(ss,:) = dataIn.data.studyRTsQuants;
    data.testRTsQuants(ss,:)  = dataIn.data.testRTsQuants;
    
    % get data
    for ff = 1:numel(itcFields)
        data.(itcFields{ff})    = cat(1,data.(itcFields{ff}),dataIn.data.(itcFields{ff}));
    end
    % get additional info
    if ss==1
        data.SR             = dataIn.data.SR;
        data.trialTime      = dataIn.data.trialTime;
    end
    
    % load pca glm data
    dataIn = load([opts.dataPath subjects{ss} '/ITC_Data/' fileName2 ]);
    for ff = 1:numel(itcGLMFields)
        data.(itcGLMFields{ff})     = cat(1,data.(itcGLMFields{ff}),dataIn.data.(itcGLMFields{ff}));
    end
    % get pca-glm additional info
    if ss==1
        data.NumComponents = dataIn.data.NumComponents;
    end

     % load MI data
    dataIn = load([opts.dataPath subjects{ss} '/MI_Data/' fileName3 ]);

    
    gitfor ff = 1:numel(miFields)
        data.(miFields{ff})     = cat(1,data.(miFields{ff}),dataIn.data.(miFields{ff}));
    end
    % get additional info
    if ss==1
        data.nPhBins = dataIn.data.nPhBins;
    end
end

data.nChans         = numel(data.ROIid);
data.nTrialSamps    = numel(data.trialTime);
data.nFreq          = nFreqs;