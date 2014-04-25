function [teststamps hcmf hmcfn hmcfnRT testRTs pairs S] = SS3_analyze_ecog(blocknum,sub, path)
%ecog behavioral data analysis
%a.gonzl
%Aug 5,2011
%mapping of responses
% sub1 button pressed   '+' ->  [2] new
%                       '4' ->  [1] old
%
% subject 2 had enough trials with confidence
% sub2 button pressed   '+' ->  [1] HC old   % switched
%                       '6' ->  [1] HC old
%                       '5' ->  [2] LC new
%                       '4' ->  [2] HC new
% sub3 button pressed   '+' ->  [2] new
%                       '4' ->  [1] old
% events codes:
% 0    -> no answer or late answer
% 1     -> HChits
% 2     -> HCcr
% 3     -> miss
% 4     -> fa
% 5     -> LChits
% 6     -> LCcr
% 7     -> no answer to old stim
% 8     -> no answer to new stim


% load data

%cd(['/biac4-wagner/biac3/wagner7/ecog/subj' num2str(subjnum) '/ecog/SS3/BehavData/'])

if nargin<3
    blocknum = 0;
end


datapath = [path sprintf('AGtestlonglist.%s.out.mat',sub)];
load(datapath)

[teststamps hcmf hmcfn hmcfnRT testRTs pairs S] = SS3_analyze_test(theData,sub,blocknum);

end

function [stamps hcmf hmcfn hmcfnRT RT pairs S] = SS3_analyze_test(data,sub,blocknum)

triallenght = numel(data.item);
numofblocks =4;
blocktrials = 1:triallenght; % all blocks

if blocknum ~=0
    blocktrials = reshape(blocktrials,[],numofblocks);
    blocktrials = blocktrials(:,blocknum);
end

% loads conditions
cond = data.cond(blocktrials);

% collapses across responses made during stim time and post stim
% time

resp = cell(numel(blocktrials),1);
RT = NaN(numel(blocktrials),1);
stamps = NaN(numel(blocktrials),1);
cnt = 1;
for trial = reshape(blocktrials,1,[])
    try
        % check for response at stim time
        if ~strcmp(data.stimresp(trial),'n')
            resp(cnt,1) = data.stimresp(trial);
            RT(cnt,1) = data.stimRT(trial);
            % check for response at post stim time
        elseif ~strcmp(data.poststimresp(trial), 'noanswer') ...
                && ~strcmp(data.poststimresp(trial), 'n');
            resp(cnt,1) = data.poststimresp(trial);
            RT(cnt,1) = data.poststimRT(trial)+1;
            % no response or other invalid response
        else
            resp(cnt,1) = {'n'};
            RT(cnt,1) = NaN;
        end
    catch ME
        display(ME)
        resp(cnt,1)=0;
    end
    if iscell(resp{cnt}) && numel(resp{cnt})==1
        resp(cnt,1) = resp{cnt};
    end
    stamps(cnt) = data.flip(trial).VBLTimestamp;
    cnt = cnt+1;
end

% pairs contains the response and condition pairs
pairs = [resp num2cell(cond)];

% mapping of button presses
if strcmp(sub,'MD')
    oldbutton = {'4','5'}; newbutton = {'+','6'};
elseif strcmp(sub,'sr') && run == 4
    oldbutton = {'L','4'}; newbutton = {'+','6'};
elseif strcmp(sub,'sr') || strcmp(sub,'RB')
    oldbutton = {'4','5'}; newbutton = {'+','6'};
elseif strcmp(sub,'DZ')
    oldbutton = {'L','5'}; newbutton = {'+','R'};
else
    oldbutton = {'+','6'}; newbutton = {'4','5'};
end

S.old = cell2num(pairs(:,2))==1; % old stimuli
S.new = cell2num(pairs(:,2))==2; % new stimuli
S.HCOresp = strcmp(oldbutton{1},pairs(:,1)); % high confidence old response
S.HCNresp = strcmp(newbutton{1},pairs(:,1)); % high confidence new response
S.LCNresp = strcmp(newbutton{2},pairs(:,1)); % low confidence new response
S.LCOresp = strcmp(oldbutton{2},pairs(:,1)); % low confidence old response
S.NoAn = strcmp('n',pairs(:,1)); % no response

S.NoAn2Old = S.NoAn & S.old;
S.NoAn2New = S.NoAn & S.new;

% hits calculation
S.HChits = S.HCOresp & S.old;
S.LChits = S.LCOresp & S.old;

% cr calculation
S.HCcr = S.HCNresp & S.new;
S.LCcr = S.LCNresp & S.new;

% trying to correct DZ
% for all of the no responses, assign to LCOresp
if strcmp(sub,'DZ')
    LCOresp = S.NoAn;
    S.modificationforDZ = 'for all of the no responses, assign to LCOresp';
    S.LChits = S.LChits| (LCOresp & S.old);
    
    S.hits = S.HChits | S.LChits;
    S.miss = ~S.hits & S.old;
    S.cr = S.HCcr | S.LCcr;
    S.fa = ~S.cr & S.new ;

else
    
    S.hits = S.HChits | S.LChits;
    S.cr = S.HCcr | S.LCcr;
    S.fa = ~S.cr & S.new & ~S.NoAn;
    S.miss = ~S.hits & S.old & ~S.NoAn;
end


% performance measures
perf_mean = mean(sum([S.hits S.cr],2));
%perf_std = std(sum([S.hits S.cr],2));

hcmf = zeros(numel(blocktrials),1);
hcmf(S.HChits) = 1;
hcmf(S.HCcr) = 2;
hcmf(S.miss) = 3;
hcmf(S.fa) = 4;
hcmf(S.LChits) = 5;
hcmf(S.LCcr) = 6;
hcmf(S.NoAn2Old) = 7;
hcmf(S.NoAn2New) = 8;


S.nHChits = sum(S.HChits);
S.nLChits = sum(S.LChits);
S.nhits = sum(S.hits);
S.nmiss = sum(S.miss);
S.nHCcr = sum(S.HCcr);
S.nLCcr = sum(S.LCcr);
S.ncr = sum(S.cr);
S.nfa = sum(S.fa);
S.noan = sum(S.NoAn);
S.nNoAn2Old = sum(S.NoAn2Old);
S.nNoAn2New = sum(S.NoAn2New);

hitsRT = mean(RT(S.hits));
S.hitRTs = RT(S.hits);

missRT = mean(RT(S.miss));
S.missRTs = RT(S.miss);

crRT = mean(RT(S.cr));
S.crRTs = RT(S.cr);

faRT = mean(RT(S.fa));
S.faRTs = RT(S.fa);

S.allRTs = RT;



% performance info
if blocknum ~=0
    fprintf('\nsubject: %s block#: %i \n', sub, blocknum)
else
    fprintf('\nsubject: %s \n', sub)
end
fprintf(['hits  ' num2str(S.nhits) '\n']);
fprintf(['miss  ' num2str(S.nmiss) '\n']);
fprintf(['corj  ' num2str(S.ncr) '\n']);
fprintf(['flsa  ' num2str(S.nfa) '\n']);
fprintf(['noan  ' num2str(S.noan) '\n']);
fprintf(['perf mean  ' num2str(perf_mean) '\n']);

[dPrime c] = calc_dPrime(S.nhits, S.nmiss, S.nfa, S.ncr);

if dPrime > 1e3
    [dPrime c] = calc_dPrime(nhits+1/triallenght, nmiss+1/triallenght, nflsa+1/triallenght, ncorj+1/triallenght);
end


fprintf(['dprime: ' num2str(dPrime) '\n']);
fprintf(['c: ' num2str(c) '\n\n']);

hmcfn = [S.nhits S.nmiss S.ncr S.nfa S.noan];
hmcfnRT = [hitsRT missRT crRT faRT 0];

fprintf('RTs (hits miss crj flsa)\n');
disp([hitsRT missRT crRT faRT])

% stamps = [data.flip.VBLTimestamp]';
return

end