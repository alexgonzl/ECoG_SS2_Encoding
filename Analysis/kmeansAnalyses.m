function out = kmeansAnalyses(data,opts)
% function that uses KMeans to find electrodes that match in their spectro-temporal patterns
% note that this is working at the T mean channel level. 
% possible alternative approaches:
%
% data must be the output of groupLPCDataMultiBand()
%	depending on analysis the input will change, however they should 
% 	have the following formatting:
% 	bands X channels X bins
%	
%	other required inputs:
% 	AnalysisBins-> this input, subselects time bins for the analysis.
%	
% required inputs in opts:
% 	chans -> vector of channels to be grouped.
%	analysis 
%		-> 'activity', clusters based on the amplitude of the signals
%		-> 'studyRT',  clusters based on the correlation of amplitude
%						and encoding reaction time
% 		-> 'testRT', 	clusters based on the correlation of 
% 					amplitude and retrieval reaction time 
% 	numClusters -> number of clusters for the kmeans analysis
% 	replicates  -> number of replications for the analysis
%
% Alex G. 
% Last updated 3/23/16
% 
%% 	
out = [];
out.opts  = opts;

% select channels based on criterion
out.chans = find(opts.chans);
out.AnalysisBins  = data.AnalysisBins;
out.nBands = data.nBands;
out.nAnalysisBins = sum(out.AnalysisBins);
switch opts.analysis
    case 'activity'
        X = permute(data.tBinChResp,[2 1 3]);        
    case 'studyRT'
        X=permute(data.dataToStudyRTsCorr,[2 1 3]);
    case 'testRT'
        X=permute(data.dataToTestRTsCorr,[2 1 3]);
end
X = X(out.chans,:,out.AnalysisBins);
X = X(:,:);

%%
out.X = X;
out.K = opts.numClusters;
[out.IDX, out.C, out.SUMD, out.D] =...
    kmeans(X, out.K,'replicates',opts.replicates,'distance','correlation');

% use original time-course for getting summary statitics (for
% display purposes)
out.Bins = data.Bins;
out.tBinCluster_tChResp = zeros(out.K, out.nBands,data.nBins);
out.pBinCluster_tChResp = zeros(out.K, out.nBands,data.nBins);
for kk = 1:out.K
    chans = out.chans(out.IDX==kk);    
    for bb = 1:out.nBands
        X =  squeeze(data.tBinChResp(bb,chans,:));
        [~,p,~,t] = ttest(X);
        out.tBinCluster_tChResp(kk,bb,:)=t.tstat;
        out.pBinCluster_tChResp(kk,bb,:)=p;
    end
end
% %% 
% % Get K-Means fot testCorrs
% X=permute(data.dataToTestRTsCorr,[2 1 3]);
% data.KMeansTest.K = 3;
% [data.KMeansTest.IDX, data.KMeansTest.C, data.KMeansTest.SUMD, data.KMeansTest.D]  = kmeans(X(:,:), data.KMeansTest.K,'replicates',100,'distance','correlation');
% data.KMeansTest.ReshapeC = zeros(data.KMeansTest.K,data.nBands,data.nBins);
% data.KMeansTest.TTests = zeros(3,data.nBands,data.nBins);
% data.KMeansTest.CorrMeans = zeros(3,data.nBands,data.nBins);
% for kk = 1:data.KMeansTest.K
%     data.KMeansTest.ReshapeC(kk,:,:) = reshape(data.KMeansTest.C(kk,:),[data.nBands,data.nBins]);
%     chans = data.KMeansTest.IDX==kk;
%     for ba= 1:data.nBands
%         temp=squeeze(data.dataToTestRTsCorr(ba,chans,:));
%         data.KMeansTest.CorrMeans(kk,ba,:) = mean(temp);
%         [~,~,~,t] = ttest(temp);        
%         data.KMeansTest.TTests(kk,ba,:)=t.tstat;
%     end 
% end
% 
% %%
% % Get K-Means fot studyCorrs
% 
% data.KMeansStudy.K = 3;
% [data.KMeansStudy.IDX, data.KMeansStudy.C, data.KMeansStudy.SUMD, data.KMeansStudy.D] ...
%     = kmeans(X(:,:), data.KMeansStudy.K,'replicates',100,'distance','correlation');
% data.KMeansStudy.ReshapeC = zeros(data.KMeansStudy.K,data.nBands,data.nBins);
% data.KMeansStudy.TTests = zeros(3,data.nBands,data.nBins);
% data.KMeansStudy.CorrMeans = zeros(3,data.nBands,data.nBins);
% for kk = 1:data.KMeansStudy.K
%     data.KMeansStudy.ReshapeC(kk,:,:) = reshape(data.KMeansStudy.C(kk,:),[data.nBands,data.nBins]);
%     chans = data.KMeansStudy.IDX==kk;
%     for ba= 1:data.nBands
%         temp=squeeze(data.dataToStudyRTsCorr(ba,chans,:));
%         data.KMeansStudy.CorrMeans(kk,ba,:) = mean(temp);
%         [~,~,~,t] = ttest(temp);
%         data.KMeansStudy.TTests(kk,ba,:)=t.tstat;
%     end 
% end