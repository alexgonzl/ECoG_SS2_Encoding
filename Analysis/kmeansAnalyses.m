function out = kmeansAnalyses(data,opts)
% function that uses KMeans to find electrods shtat match thei spectro-temporal patterns.
% note that this is working at the T mean channel level. 
% possible alternative approaches:
% 1) work at the mean channel activity level 
% 2) pre-select channels based on some criteria (like overall channel activity through trials/time/frequency)
% 3) do independent groupings of channels per analysis; or carry the same groupings into the rest of the analyses
% 4) somehow sample trials and create a distribution of channel groupings. 
%
% data must be the output of groupLPC data

%% 
% select channels based on criterion
channels = sum(data.BinChRespScore>opts.opts.TvalChanThr)>1;

switch opts.analysis
    case 'activity'
        X = permute(data.BinERPsAvg,[2 1 3]);        
    case 'studyRT'
        X=permute(data.dataToStudyRTsCorr,[2 1 3]);
    case 'testRT'
        X=permute(data.dataToTestRTsCorr,[2 1 3]);

end
X = X(channels,:,:);

%%
% Get K-Means for activity

data.KMeansERA.K = 3;
[data.KMeansERA.IDX, data.KMeansERA.C, data.KMeansERA.SUMD, data.KMeansERA.D]  = kmeans(X(:,:), data.KMeansERA.K,'replicates',100,'distance','correlation');
data.KMeansERA.ReshapeC = zeros(data.KMeansERA.K,data.nBands,data.nBins);
data.KMeansERA.TTests = zeros(3,data.nBands,data.nBins);
for kk = 1:data.KMeansERA.K
    data.KMeansERA.ReshapeC(kk,:,:) = reshape(data.KMeansERA.C(kk,:),[data.nBands,data.nBins]);        
    chans = data.KMeansERA.IDX==kk;
    for ba= 1:data.nBands
        temp=squeeze(data.BinERPsAvg(ba,chans,:));
        [~,~,~,t] = ttest(temp);
        data.KMeansERA.TTests(kk,ba,:)=t.tstat;
    end 
end

%% 
% Get K-Means fot testCorrs
X=permute(data.dataToTestRTsCorr,[2 1 3]);
data.KMeansTest.K = 3;
[data.KMeansTest.IDX, data.KMeansTest.C, data.KMeansTest.SUMD, data.KMeansTest.D]  = kmeans(X(:,:), data.KMeansTest.K,'replicates',100,'distance','correlation');
data.KMeansTest.ReshapeC = zeros(data.KMeansTest.K,data.nBands,data.nBins);
data.KMeansTest.TTests = zeros(3,data.nBands,data.nBins);
data.KMeansTest.CorrMeans = zeros(3,data.nBands,data.nBins);
for kk = 1:data.KMeansTest.K
    data.KMeansTest.ReshapeC(kk,:,:) = reshape(data.KMeansTest.C(kk,:),[data.nBands,data.nBins]);
    chans = data.KMeansTest.IDX==kk;
    for ba= 1:data.nBands
        temp=squeeze(data.dataToTestRTsCorr(ba,chans,:));
        data.KMeansTest.CorrMeans(kk,ba,:) = mean(temp);
        [~,~,~,t] = ttest(temp);        
        data.KMeansTest.TTests(kk,ba,:)=t.tstat;
    end 
end

%%
% Get K-Means fot studyCorrs

data.KMeansStudy.K = 3;
[data.KMeansStudy.IDX, data.KMeansStudy.C, data.KMeansStudy.SUMD, data.KMeansStudy.D] ...
    = kmeans(X(:,:), data.KMeansStudy.K,'replicates',100,'distance','correlation');
data.KMeansStudy.ReshapeC = zeros(data.KMeansStudy.K,data.nBands,data.nBins);
data.KMeansStudy.TTests = zeros(3,data.nBands,data.nBins);
data.KMeansStudy.CorrMeans = zeros(3,data.nBands,data.nBins);
for kk = 1:data.KMeansStudy.K
    data.KMeansStudy.ReshapeC(kk,:,:) = reshape(data.KMeansStudy.C(kk,:),[data.nBands,data.nBins]);
    chans = data.KMeansStudy.IDX==kk;
    for ba= 1:data.nBands
        temp=squeeze(data.dataToStudyRTsCorr(ba,chans,:));
        data.KMeansStudy.CorrMeans(kk,ba,:) = mean(temp);
        [~,~,~,t] = ttest(temp);
        data.KMeansStudy.TTests(kk,ba,:)=t.tstat;
    end 
end