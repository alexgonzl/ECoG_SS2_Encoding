function S = SS2_BehavAnalysis(S,run,study,test)

nItems   = numel(study.item);
% get indices of study items in the test data
[~,b]   = ismember(study.item,test.item);

tKeys = S.TestKeys(run,:);
sKeys = S.StudyKeys(run,:);

studyCond   = study.cond;
testCond    = test.cond(b);

studyRTs    = nan(nItems,1);
testRTs     = nan(nItems,1);
for ii = 1:nItems
    if study.respRT(ii)>0
        studyRTs(ii)   = study.respRT(ii);
    end
    if test.stimRT(b(ii)) > 0
        testRTs(ii)     = test.stimRT(b(ii));
    elseif test.poststimRT(b(ii)) > 0
        testRTs(ii)     = test.poststimRT(b(ii))+1;
    end       
end
S.studyRTs  = [S.studyRTs; studyRTs];
S.testRTs   = [S.testRTs ; testRTs];

tResp = nan(nItems,1);
for  jj = 1:numel(tKeys)
    [~,idx]     = ismember(test.stimresp(b),tKeys{jj});
    tResp(idx==1)  = jj;
    [~,idx]     = ismember(test.poststimresp(b),tKeys{jj});
    tResp(idx==1)  = jj;
end

assert(mean(isnan(testRTs)==isnan(tResp)), 'respones don''t match the recorded RTs') 

S.old = cell2num(pairs(:,2))==1; % old stimuli
S.new = cell2num(pairs(:,2))==2; % new stimuli
S.HCOresp = strcmp(oldbutton{1},pairs(:,1)); % high confidence old response
S.LCOresp = strcmp(oldbutton{2},pairs(:,1)); % low confidence old response
S.HCNresp = strcmp(newbutton{1},pairs(:,1)); % high confidence new response
S.LCNresp = strcmp(newbutton{2},pairs(:,1)); % low confidence new response
S.NoAn = strcmp('n',pairs(:,1)); % no response

% hits calculation
S.HChits = S.HCOresp & S.old;
S.LChits = S.LCOresp & S.old;
S.hits = S.HChits | S.LChits;
S.miss = ~S.hits & S.old & ~S.NoAn;
S.NoAn2Old = S.NoAn & S.old;

% cr calculation
S.HCcr = S.HCNresp & S.new;
S.LCcr = S.LCNresp & S.new;
S.cr = S.HCcr | S.LCcr;
S.fa = ~S.cr & S.new & ~S.NoAn;
S.NoAn2New = S.NoAn & S.new;