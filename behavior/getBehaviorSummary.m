% SS2e Analyze Behavior

behavPerf = [];
behavPerf.subjects      = {'SRb','MD','LK','NC'};
behavPerf.dataPath      = '/Volumes/ECoG_SS2/SS2/SS2e/Results/';
nItemsPerBlock = 20;

preAlloc              = {'condID','sRespID','tRespID','blockID','subjID','studyRTs','testRTs'};
for aa = 1:numel(preAlloc)
    behavPerf.(preAlloc{aa}) = [];
end

nSubjs = numel(behavPerf.subjects);
for ss = 1:nSubjs
    behavDataPath = [behavPerf.dataPath behavPerf.subjects{ss} '/BehavResults/'];
    load([behavDataPath 'behavResults.mat']);
    
    
    conds   = unique(S.cond);
    nConds  = numel(conds);
    for jj = 1:nConds
        behavPerf.Subj_studyRT_EncPerf(ss,jj)  =  getVectorStats(S.studyRTs(S.cond==conds(jj)));
        behavPerf.Subj_testRT_EncPerf(ss,jj)   =  getVectorStats(S.testRTs(S.cond==conds(jj)));
    end
    
    % works for two conditions
    behavPerf.SubjStudyPerf(ss)         = getBinaryPerf(S.cond,S.sResp,conds(1),conds(2));
    behavPerf.Subj_studyRTCondDiff(ss)  = getTwoSampStats(S.studyRTs(S.cond==conds(1)),S.studyRTs(S.cond==conds(2)));
    behavPerf.Subj_testRTCondDiff(ss)   = getTwoSampStats(S.testRTs(S.cond==conds(1)),S.testRTs(S.cond==conds(2)));
    
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
