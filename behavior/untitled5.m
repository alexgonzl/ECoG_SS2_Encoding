
% Analysis of behavioral data for SS2 script
% Apr 2014
% A. Gonzl

% dependencies:
% ----

S =[];
S.expt  = 'SS2e';
S.subjName   = 'SRb';
S = SS2e_subjInfo(S);
S.behavDataPath = ['/biac4/wagner/biac3/wagner7/ecog/subj' S.subjNum '/ecog/SS2/BehavData/'];

StudyConds = {'Abs','Conc','CorrectAbs','InCorrectAbs', 'CorrectConc','InCorrectConc','noAnswer'};
for i=1:numel(StudyConds)
    S.sConds.(StudyConds{i}) = [];
end

otherFields = {'items','cond','sResp','tResp','subRem','subForg','studyRTs','testRTs'};
for i=1:numel(otherFields)
    S.(otherFields{i}) = [];
end

for run = S.run_nums
    study   = load( [S.behavDataPath 'study_' num2str(run) '.' S.subjName '.out.mat']);
    study   = study.theData;
    test    = load( [S.behavDataPath 'test_' num2str(run) '.' S.subjName '.out.mat']);
    test    = test.theData;
   
    S = SS2e_CodeTrials(S,run,study,test);
end
% 
% for run = 1: S.nruns
%     run_num = S.run_nums(run);
%     
%     [teststamps hcmf hmcfn hmcfnRT testRTs pairs{run} S.S{run}] = SS2_analyze_ecog(run_num,S.sub,path);
%     
%     
%     for i=1:numel(S.conds)
%         S.(S.conds{i}) = vertcat(S.(S.conds{i}), ...
%             S.S{run}.(S.conds{i}));
%     end
%     
%     %S.codes = vertcat(S.codes, hcmf);
%     
% end
% 
% for i=1:numel(S.conds)
%     S.nConds.(S.conds{i}) = sum(S.(S.conds{i}));
% end
% 
% [S.dPrime S.c] = calc_dPrime(sum(S.nhits), sum(S.nmiss), sum(S.nfa), sum(S.ncr));
% S.hitrate = sum(S.nhits) / (sum(S.nhits)+ sum(S.nmiss));
% S.farate = sum(S.nfa) / (sum(S.nfa) + sum(S.ncr));
% 
% S.hitRTStats = [];
% S.hitRTStats = calc_RTstats(S.hitRTs,S.hitRTStats);
% 
% S.crRTStats = [];
% S.crRTStats = calc_RTstats(S.crRTs,S.crRTStats);
% 
% S.missRTStats = [];
% S.missRTStats = calc_RTstats(S.missRTs,S.missRTStats);
% 
% S.faRTStats = [];
% S.faRTStats = calc_RTstats(S.faRTs,S.faRTStats);
% 
% % ttest
% [H,P,CI] = ttest2(S.hitRTs,S.crRTs);
% S.two_samp_ttest_P_hcr = P;
% S.two_samp_ttest_CI_hcr = CI;
% 
% [H,P,CI] = ttest2(S.hitRTs,S.faRTs);
% S.two_samp_ttest_P_hfa = P;
% S.two_samp_ttest_CI_hfa = CI;
% 
% save([path '/S.mat'],'S')
% path2= '/Users/alexgonzalez/Documents/ECOG/Results/BehavData/';
% save([path2 'Subj' S.subnum 'S.mat'],'S')
% 
% fprintf('\n\nTotal Performance for subject: %s \n', S.sub)
% fprintf(['Hit Rate : ' num2str(S.hitrate) '\n']);
% fprintf(['FA Rate  : ' num2str(S.farate) '\n']);
% fprintf(['dprime   : ' num2str(S.dPrime) '\n']);
% fprintf(['c        : ' num2str(S.c) '\n\n']);
% 
% fprintf('RTs (hits miss crj flsa)\n');
% disp([S.hitRTStats.RTmedian S.missRTStats.RTmedian ...
%     S.crRTStats.RTmedian S.faRTStats.RTmedian])
% fprintf('P value for hitRTs vs crRTs\n')
% fprintf(['P         : ' num2str(S.two_samp_ttest_P_hcr) '\n\n'])
% 
% fprintf('P value for hitRTs vs faRTs\n')
% fprintf(['P         : ' num2str(S.two_samp_ttest_P_hfa) '\n\n'])
% 
% 
% 
% 
