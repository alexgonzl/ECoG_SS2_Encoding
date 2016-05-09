% Analysis of behavioral data for SS2 script
% also stores markers for pre-proecessing
% Apr 2014
% A. Gonzl

% dependencies:
% SS2e_subjInfo
% SS2e_CodeTrials
% findEventMarkers
subjects      = {'16b','17b','18','19','24','28','29','30'};

for ss = 1:numel(subjects);
    S =[];
    S.expt  = 'SS2e';
    S.subjNum   = subjects{ss};
    S = SS2e_subjInfo(S,S.subjNum);
    
    S.dataPath      = ['/Volumes/ECoG_SS2/SS2/data/' S.subjName '/' ];
    S.behavDataPath = [S.dataPath 'BehavData/'];
    S.RawDataPath   = [S.dataPath 'RawData/'];
    
    S.behavResPath  = ['/Volumes/ECoG_SS2/SS2/SS2e/Results/' S.subjNum '/BehavResults/' ];
    S.behavResPath2 = ['~/Google Drive/Research/ECoG_SS2e/Behavior/'  S.subjNum '/'] ;
    
    StudyConds = {'Abs','Conc','CorrectAbs','InCorrectAbs', 'CorrectConc','InCorrectConc','noAnswer'};
    for i=1:numel(StudyConds)
        S.sConds.(StudyConds{i}) = [];
    end
    
    otherFields = {'items','cond','sResp','tResp','subRem','subForg','studyRTs','testRTs','studyIDatTest'};
    for i=1:numel(otherFields)
        S.(otherFields{i}) = [];
    end
    
    cnt = 1;
    for run = S.run_nums
        fprintf( '\nProcessing Subject %s, Block %s, Run# %d \n',S.subjName, S.blocklist{cnt}, run);
        if ~strcmp(S.subjNum,'19')
            study   = load( [S.behavDataPath 'study_' num2str(run) '.' S.subjName '.out.mat']);
            study   = study.theData;
            test    = load( [S.behavDataPath 'test_' num2str(run) '.' S.subjName '.out.mat']);
            test    = test.theData;
        else
            study   = load( [S.behavDataPath 'AGstudylonglist.' S.subjName '.out.mat']);
            study   = study.theData;
            test    = load( [S.behavDataPath 'AGtestlonglist.' S.subjName '.out.mat']);
            test    = test.theData;
        end
        
        S = SS2e_CodeTrials(S,cnt,study,test);
        
        % get event time stamps
        load([S.RawDataPath S.blocklist{cnt} '/Pdio' S.blocklist{cnt} '_02.mat']);
        S.eventTimeStamps{cnt} = findEventMarkers(anlg, [study.stimFlip.VBLTimestamp]);
        cnt = cnt+1;
    end
    
    if ~exist(S.behavResPath,'dir'), mkdir(S.behavResPath),end
    save([S.behavResPath 'behavResults'], 'S')
    save([S.behavResPath2 'behavResults'], 'S')
end
