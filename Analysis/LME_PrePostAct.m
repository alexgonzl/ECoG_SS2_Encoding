

addpath lib/fdr_bh/
% re-organize data for model analyses.


% Study
Y1=data.StudyRTs_PrePostActModel_Ts;

% bands x channel x pre/post/inter
[nBands,nChans,nConds] = size(Y1);
N = numel(Y1);
study = Y1(:); 

% Test
test=data.TestRTs_PrePostActModel_Ts(:);

% set up of independent variables to match data:
% frequency bands
bands = repmat((1:nBands)',[nChans*nConds,1]);
bands_IDs = {'delta','theta','alpha','beta','lgam','hgam'};
bands_s = cell(N,1);
for ii = 1:6
   bands_s(bands==ii)= bands_IDs(ii);
end

% conditions pre=1,post=2, inter=3
conds = [ones(nBands*nChans,1);2*ones(nBands*nChans,1);3*ones(nBands*nChans,1)];

conds_s = cell(N,1);
conds_s(conds==1)={'pre'};
conds_s(conds==2)={'post'};
conds_s(conds==3)={'post-pre'};

% channels
temp=repmat(1:nChans,[nBands,1]);
chans=repmat(temp(:),[nConds,1]);

% channels by subject ID 
nSubjs=8;
subjChans = zeros(N,1);
for ii=1:nSubjs
    chanIDs = find(data.subjChans==ii);
    subjChans(ismember(chans,chanIDs))=ii;    
end

% channels by hemisphere
hemChans = zeros(N,1);
for ii=1:2
    chanIDs = find(data.hemChan==ii);
    hemChans(ismember(chans,chanIDs))=ii;
end

% region: ips, spl, ag
rois = zeros(N,1);
rois_s = cell(N,1);
rois_ids = {'IPS','SPL','AG'};
for ii=1:3
    chanIDs = find(data.ROIid==ii);
    rois(ismember(chans,chanIDs))=ii;
    rois_s(ismember(chans,chanIDs))=rois_ids(ii);
end

tbl2=table(study,test,bands_s, conds_s, ...
    rois_s, categorical(subjChans), categorical(hemChans),...
    'variablenames',{'study','test','bands','conds','rois','subj','hem'});

X=grpstats(tbl2,{'bands','conds'},{'mean','sem'},'datavars',{'study','test'});

testRT_Tvals_CondsBand = zeros(nConds,nBands,2);
studyRT_Tvals_CondsBand = zeros(nConds,nBands,2);

testRT_Tvals_CondsBandROIs = zeros(nConds,nBands,3,2);
studyRT_Tvals_CondsBandROIs = zeros(nConds,nBands,3,2);

for ii = 1:nConds
    for jj = 1:nBands
        [~,p,~,t] = ttest(test(conds==ii & bands == jj));
        testRT_Tvals_CondsBand(ii,jj,1) = t.tstat;
        testRT_Tvals_CondsBand(ii,jj,2) = p;
        
        [~,p,~,t] = ttest(study(conds==ii & bands == jj));
        studyRT_Tvals_CondsBand(ii,jj,1) = t.tstat;
        studyRT_Tvals_CondsBand(ii,jj,2) = p;
        for kk=1:3
            [~,p,~,t] = ttest(test(conds==ii & bands == jj & rois==kk));
            testRT_Tvals_CondsBandROIs(ii,jj,kk,1) = t.tstat;
            testRT_Tvals_CondsBandROIs(ii,jj,kk,2) = p;
            
            [~,p,~,t] = ttest(study(conds==ii & bands == jj & rois==kk));
            studyRT_Tvals_CondsBandROIs(ii,jj,kk,1) = t.tstat;  
            studyRT_Tvals_CondsBandROIs(ii,jj,kk,2) = p;          
        end
    end
end

%% study model
disp('Full model for study data')
lmeFull = fitlme(tbl2,'study~1+bands*conds*rois+hem + (1|subj)');
disp(anova(lmeFull))
% compute if reduce model is significant
lmeInter = fitlme(tbl2,'study~1+bands*conds*rois-bands:conds:rois +hem + (1|subj)');
compare(lmeInter,lmeFull)
%% test model
disp('Full model for test data')
lmeFull = fitlme(tbl2,'test~1+bands*conds*rois+hem + (1|subj)');
disp(anova(lmeFull))
% compute if reduce model is significant
lmeInter = fitlme(tbl2,'test~1+bands*conds*rois-bands:conds:rois +hem + (1|subj)');
compare(lmeInter,lmeFull)

