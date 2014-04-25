function S = SS2e_subjInfo(S)
% function returns relevant information for the encodign data for each subject

S.StudyResponses = {'abs','conc'};
S.TestResposnes  = {'HCOresp','LCOresp','LCNresp','HCNresp'};
switch S.subjName
    case {'16b','SRb'}   
    S.subjNum = '16b';
    S.nruns = 4;
    S.run_nums = 1:4;
    S.blocklist = {'SRb-25', 'SRb-27', 'SRb-29', 'SRb-31'};    
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
end

% elseif strcmp(subj,'RHb')
%     S.subnum = '17b';
%     S.nruns = 2;
%     S.run_nums = 1:2;
% elseif strcmp(subj,'MD')
%     subnum = '18';
%     S.nruns=2;
%     S.run_nums = 1:2;
% elseif strcmp(subj,'TC')
%     subnum = '22';
%     S.nruns=3;
%     S.run_nums = 1:3;
% elseif strcmp(subj,'LK')
%     subnum = '24';
%     S.nruns=4;
%     S.run_nums = 2:5;
% elseif strcmp(subj,'KS')
%     subnum = '27';
%     S.nruns=4;
%     S.run_nums = 1:4;
% elseif strcmp(subj,'NC')
%     subnum = '28';
%     S.nruns=4;
%     S.run_nums = 1:4;
% elseif strcmp(subj,'JT')
%     subnum = '29';
%     S.nruns=4;
%     S.run_nums = 1:4;
% end