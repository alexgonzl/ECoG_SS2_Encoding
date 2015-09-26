function S = SS2e_CodeTrials(S,run,study,test)
% coding and storing of responses for the behavioral data 

nItems   = numel(study.item);
% get indices of study items in the test data
[~,b]   = ismember(study.item,test.item);

% store study IDs at test, and test IDs at study.
S.studyIDatTest = [S.studyIDatTest; b+(run-1)*numel(test.item)];
% update all the test data to only the ones that were at encoding, i.e. old
% items
testFields = {'item','cond','stimRT','poststimRT','stimresp','poststimresp'};

for ii = 1:numel(testFields)
    test.(testFields{ii}) = test.(testFields{ii})(b);
end

tKeys = S.TestKeys(run,:);
sKeys = S.StudyKeys(run,:);

studyCond   = study.cond;

S.items = [S.items;study.item];
S.cond  = [S.cond;studyCond];

% store RTs
studyRTs    = nan(nItems,1);
testRTs     = nan(nItems,1);
for ii = 1:nItems
    if study.respRT(ii)>0
        studyRTs(ii)   = study.respRT(ii);
    elseif study.ITIrespRT(ii)>0
        studyRTs(ii)   = study.ITIrespRT(ii)+2.5; % 2.5s= stimTime
    end
    if test.stimRT(ii) > 0
        testRTs(ii)     = test.stimRT(ii);
    elseif test.poststimRT(ii) > 0
        testRTs(ii)     = test.poststimRT(ii)+1;% 1 = stimTime
    end       
end
S.studyRTs  = [S.studyRTs; studyRTs];
S.testRTs   = [S.testRTs ; testRTs];

% store responses
sResp   = codeCondition(study.resp,sKeys);
temp    = codeCondition(study.ITIresp,sKeys);
sResp(~isnan(temp)) = temp(~isnan(temp));

tResp   = codeCondition(test.stimresp,tKeys);
temp    = codeCondition(test.poststimresp,tKeys);
tResp(~isnan(temp)) = temp(~isnan(temp));

%assert(mean(isnan(studyRTs)==isnan(sResp))==1, 'respones don''t match the recorded RTs') 
%assert(mean(isnan(testRTs)==isnan(tResp))==1,  'respones don''t match the recorded RTs') 

% code responses
S.sResp     = [S.sResp  ; sResp];
S.tResp     = [S.tResp  ; tResp];

S.subRem    = [S.subRem     ; (tResp==1 | tResp==2)]; % responded old (hit)
S.subForg   = [S.subForg    ; (tResp==3 | tResp==4)]; % responded as new (miss)

S.sConds.Abs            = [S.sConds.Abs             ; (studyCond==1)];
S.sConds.CorrectAbs     = [S.sConds.CorrectAbs      ; (sResp==1)&(studyCond==1)];
S.sConds.InCorrectAbs   = [S.sConds.InCorrectAbs    ; (sResp~=1)&(studyCond==1)];

S.sConds.Conc           = [S.sConds.Conc            ; (studyCond==2)];
S.sConds.CorrectConc    = [S.sConds.CorrectConc     ; (sResp==2)&(studyCond==2)];
S.sConds.InCorrectConc  = [S.sConds.InCorrectConc   ; (sResp~=2)&(studyCond==2)];

S.sConds.noAnswer       = [S.sConds.noAnswer        ; (sResp~=1)&(sResp~=2)];

end

function out = codeCondition(resp,keys)

nItems  = numel(resp);
out     = nan(nItems,1);
cellIdx = find(cellfun('isclass',resp,'cell'));
nCells  = numel(cellIdx);
goodIdx = setdiff(1:nItems,cellIdx);

resp1   = resp(goodIdx);
resp2   = resp(cellIdx);
for jj = 1:numel(keys)
    [~,Idx] = ismember(resp1,keys{jj});
    out(goodIdx(Idx==1))  = jj;
    
    for kk = 1:nCells
        if strcmp(resp2{kk},keys{jj})
            out(cellIdx(kk))=jj;
        end
    end
end
end

% 
% S.old = cell2num(pairs(:,2))==1; % old stimuli
% S.new = cell2num(pairs(:,2))==2; % new stimuli
% S.HCOresp = strcmp(oldbutton{1},pairs(:,1)); % high confidence old response
% S.LCOresp = strcmp(oldbutton{2},pairs(:,1)); % low confidence old response
% S.HCNresp = strcmp(newbutton{1},pairs(:,1)); % high confidence new response
% S.LCNresp = strcmp(newbutton{2},pairs(:,1)); % low confidence new response
% S.NoAn = strcmp('n',pairs(:,1)); % no response
% 
% % hits calculation
% S.HChits = S.HCOresp & S.old;
% S.LChits = S.LCOresp & S.old;
% S.hits = S.HChits | S.LChits;
% S.miss = ~S.hits & S.old & ~S.NoAn;
% S.NoAn2Old = S.NoAn & S.old;
% 
% % cr calculation
% S.HCcr = S.HCNresp & S.new;
% S.LCcr = S.LCNresp & S.new;
% S.cr = S.HCcr | S.LCcr;
% S.fa = ~S.cr & S.new & ~S.NoAn;
% S.NoAn2New = S.NoAn & S.new;