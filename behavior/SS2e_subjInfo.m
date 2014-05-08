function S = SS2e_subjInfo(S,subj)
% function returns relevant information for the encodign data for each subject

S.StudyResponses = {'abs','conc'};
S.TestResposnes  = {'HCOresp','LCOresp','LCNresp','HCNresp'};
switch subj
    case {'16b' , 'SRb'}
        S.subjNum   = '16b';
        S.subjName  = 'SRb';
        S.nruns     = 4;
        S.run_nums  = 1:4;
        S.blocklist = {'SRb-25', 'SRb-27', 'SRb-29', 'SRb-31'};
        S.StudyKeys = cell(max(S.run_nums),numel(S.StudyResponses));
        S.TestKeys = cell(max(S.run_nums),numel(S.TestResposnes));
        for rr = 1:S.nruns
            if (rr == 4)
                S.StudyKeys(rr,:) = {'R','L'};
                S.TestKeys (rr,:) = {'4','L','6','+'};
            else
                S.StudyKeys(rr,:) = {'5','4'};
                S.TestKeys (rr,:) = {'4','5','6','+'};
            end
        end
    case {'17b' , 'RHb'}
        S.subjNum = '17b';
        S.nruns = 3;
        S.run_nums = 1:3;
        S.blocklist = {'RHb0211-03', 'RHb0211-08', 'RHb0211-10'};
        S.StudyKeys = cell(max(S.run_nums),numel(S.StudyResponses));
        S.TestKeys = cell(max(S.run_nums),numel(S.TestResposnes));
        for rr = S.run_nums
            if (rr == 4)
                S.StudyKeys(rr,:) = {'R','L'};
                S.TestKeys (rr,:) = {'4','L','6','+'};
            else
                S.StudyKeys(rr,:) = {'5','4'};
                S.TestKeys (rr,:) = {'4','5','6','+'};
            end
        end
    case {'18'  , 'MD'}
        S.subjNum   = '18';
        S.subjName  = 'MD';
        S.nruns = 2;
        S.run_nums = 1:2;
        S.blocklist = {'MD0311-13', 'MD0311-16'};
        S.StudyKeys = cell(max(S.run_nums),numel(S.StudyResponses));
        S.TestKeys = cell(max(S.run_nums),numel(S.TestResposnes));
        for rr = 1:S.nruns
            S.StudyKeys(rr,:) = {'4','5'};
            S.TestKeys (rr,:) = {'4','5','6','+'};
        end
    case {'19'  , 'RB'}
        % SS3...
    case {'24'  , 'LK'}
        S.subjNum = '24';
        S.subjName = 'LK';
        S.nruns = 4;
        S.run_nums = 2:5;
        S.blocklist = {'LK_13', 'LK_15', 'LK_17', 'LK_19'};
        S.StudyKeys = cell(max(S.run_nums),numel(S.StudyResponses));
        S.TestKeys = cell(max(S.run_nums),numel(S.TestResposnes));
        for rr = 1:S.nruns
            S.StudyKeys(rr,:) = {'1','2'};
            S.TestKeys (rr,:) = {'+','6','5','4'};
        end
    case {'28'  , 'NC'}
        S.subjNum   = '28';
        S.subjName  = 'NC';
        S.nruns = 4;
        S.run_nums = 1:4;
        S.blocklist = {'NC_10', 'NC_12', 'NC_14', 'NC_16'};
        S.StudyKeys = cell(max(S.run_nums),numel(S.StudyResponses));
        S.TestKeys = cell(max(S.run_nums),numel(S.TestResposnes));
        for rr = 1:S.nruns
            S.StudyKeys(rr,:) = {'2','1'};
            S.TestKeys (rr,:) = {'4','5','6','+'};
        end
   
    case {'29'  , 'JT'}
    otherwise
        error('subject not available for SS2e')
end