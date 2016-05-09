% SS2e Analyze Behavior
addpath ~/Documents/ECoG_SS2_Encoding/lib/
savePath  = '~/Google Drive/Research/ECoG_SS2e/Behavior/summary/';

% output structure
behavPerf = [];
behavPerf.subjects      = {'16b','17b','18','19','24','28','29','30'};
behavPerf.dataPath      = '~/Google Drive/Research/ECoG_SS2e/Behavior/';

preAlloc              = {'condID','sRespID','tRespID','blockID','subjID','studyRTs','testRTs'};
for aa = 1:numel(preAlloc)
    behavPerf.(preAlloc{aa}) = [];
end

nSubjs = numel(behavPerf.subjects);
for ss = 1:nSubjs
    if strcmp(behavPerf.subjects{ss},'19')
        nItemsPerBlock = 160;
    else
        nItemsPerBlock = 20;
    end
    behavDataPath = [behavPerf.dataPath behavPerf.subjects{ss}];
    load([behavDataPath '/behavResults.mat']);
    
    
    conds   = unique(S.cond);
    nConds  = numel(conds);
    for jj = 1:nConds
        behavPerf.Subj_studyRT_EncPerf(ss,jj)  =  getVectorStats(S.studyRTs(S.cond==conds(jj)));
        behavPerf.Subj_testRT_EncPerf(ss,jj)   =  getVectorStats(S.testRTs(S.cond==conds(jj)));
    end
    [behavPerf.studytestRTsCorr(ss), behavPerf.studytestRTsPval(ss)] =  corr(S.studyRTs,S.testRTs,'rows','complete');
    
    % works for two conditions
    behavPerf.SubjStudyPerf(ss)         = getBinaryPerf(S.cond,S.sResp,conds(1),conds(2));
    behavPerf.Subj_studyRTCondDiff{ss}  = getTwoSampStats(S.studyRTs(S.cond==conds(1)),S.studyRTs(S.cond==conds(2)));
    behavPerf.Subj_testRTCondDiff{ss}   = getTwoSampStats(S.testRTs(S.cond==conds(1)),S.testRTs(S.cond==conds(2)));
    
    % store across subjects for all items
    N                   = numel(S.cond);
    nBlocks             = S.nruns;
    
    behavPerf.condID    = [behavPerf.condID     ; S.cond];
    behavPerf.sRespID   = [behavPerf.sRespID    ; S.sResp];
    behavPerf.tRespID   = [behavPerf.tRespID    ; S.tResp];
    behavPerf.blockID   = [behavPerf.blockID    ; reshape(ones(nItemsPerBlock,nBlocks)*diag(1:nBlocks),N,[])];
    behavPerf.subjID    = [behavPerf.subjID     ; ss*ones(N,1)];
    
    behavPerf.studyRTs  = [behavPerf.studyRTs   ; S.studyRTs];
    behavPerf.testRTs   = [behavPerf.testRTs     ; S.testRTs];
    
end

% works for two conditions
behavPerf.StudyPerf             = getBinaryPerf(behavPerf.condID,behavPerf.sRespID,conds(1),conds(2));
behavPerf.studyRTCondDiff       = getTwoSampStats(behavPerf.studyRTs(behavPerf.condID==conds(1)),behavPerf.studyRTs(behavPerf.condID==conds(2)));
behavPerf.testRTCondDiff        = getTwoSampStats(behavPerf.testRTs(behavPerf.condID==conds(1)),behavPerf.testRTs(behavPerf.condID==conds(2)));

save([savePath 'behavSummary.mat'],'behavPerf')
