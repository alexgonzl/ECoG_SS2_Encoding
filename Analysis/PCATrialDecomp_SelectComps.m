function out = PCATrialDecomp_SelectComps(data,opts)
% this analysis function finds the principal components for each channel
% that explains the most variance across channels.
%
% it also uses these components to
%
% opts:
% 	nComps -> number of components
%
% data input must be from the groupLPCDataMultiband.

% fix random seed for stability of clusters.
rng(1);
out                 = [];
out.nSubjs 			= size(data.BinERPs,1);
out.nComps 			= opts.nComps;
out.nChans 			= data.nChans;
out.nFeat           = data.nBands*sum(data.AnalysisBins);

% Pre allocation for outputs from PCA
out.Comps           = cell(out.nSubjs,1);
out.Projections     = zeros(out.nChans,out.nFeat,out.nComps);
out.VarExp          = zeros(out.nChans,out.nComps);

% Pre allocation for correlations
out.CorrStudyRTs    = zeros(out.nChans,out.nComps);
out.CorrTestRTs     = zeros(out.nChans,out.nComps);

% Pre allocation for GLMs
out.StudySelectedComps      = cell(out.nChans,1);
out.StudyGLMs       		= cell(out.nChans,1);
out.StudyGLMSChanCompTVal 	= nan(out.nChans,out.nComps);
out.StudyGLMsChanR2         = zeros(out.nChans,1);
out.StudyGLMsChanAR2         = zeros(out.nChans,1);

out.TestSelectedComps       = cell(out.nChans,1);
out.TestGLMs        		= cell(out.nChans,1);
out.TestGLMSChanCompTVal 	= nan(out.nChans,out.nComps);
out.TestGLMsChanR2          = zeros(out.nChans,1);
out.TestGLMsChanAR2         = zeros(out.nChans,1);

for ss=1:out.nSubjs
    subjChans = find(data.subjChans==ss);
    % re-order to channels, trials , bands , time bins
    x = permute(data.BinERPs{ss}(:,:,:,data.AnalysisBins),[2 3 1 4]);
    % concatenate bands and time bins
    x = x(:,:,:);
    [nSubjChan,nSubjTrials,~] = size(x); % d1 -> n subj channels, d2 -> # of trials, d3 -> bands*timebins
    rts1 = -log10(data.studyRTs{ss}(data.trials{ss}));
    rts2 = -log10(data.testRTs{ss}(data.trials{ss}));
    
    out.Comps{ss} = zeros(nSubjChan,nSubjTrials,out.nComps);
    
    
    % for subject channels
    for c=1:nSubjChan
        % get the trial by band*bin matrix for each channel.
        xx = squeeze(x(c,:,:));
        
        % Robust PCA to denoise; 
        A = inexact_alm_rpca(xx,0.2);
        % PCA
        [C,S,~,~,E]= pca(A','NumComponents',out.nComps);
        out.Comps{ss}(c,:,:) = C;
        out.Projections(subjChans(c),:,:) = S;
        out.VarExp(subjChans(c),:) = E(1:out.nComps);
        
        % Correlations
        r1 = corr(C,rts1,'type','spearman');
        r2 = corr(C,rts2,'type','spearman');
        out.CorrStudyRTs(subjChans(c),:) = r1;
        out.CorrTestRTs(subjChans(c),:)  = r2;
        
        CompsIdx1 = find(abs(r1)>=opts.rThr);
        out.StudySelectedComps{subjChans(c)} = CompsIdx1;
        CompsIdx2 = find(abs(r2)>=opts.rThr);
        out.TestSelectedComps{subjChans(c)} = CompsIdx2;
        
        % GLMSs
        out.StudyGLMs{subjChans(c)} = fitglm(C(:,CompsIdx1),rts1);
        out.StudyGLMSChanCompTVal(subjChans(c),CompsIdx1) ...
            = out.StudyGLMs{subjChans(c)}.Coefficients.tStat(2:end);
        out.StudyGLMsChanR2(subjChans(c))...
            = out.StudyGLMs{subjChans(c)}.Rsquared.Ordinary;
        out.StudyGLMsChanAR2(subjChans(c)) = ...
            out.StudyGLMs{subjChans(c)}.Rsquared.Adjusted;
        
        out.TestGLMs{subjChans(c)} = fitglm(C(:,CompsIdx2),rts2);
        out.TestGLMSChanCompTVal(subjChans(c),CompsIdx2) ...
            = out.TestGLMs{subjChans(c)}.Coefficients.tStat(2:end);
        out.TestGLMsChanR2(subjChans(c)) ...
            = out.TestGLMs{subjChans(c)}.Rsquared.Ordinary;
        out.TestGLMsChanAR2(subjChans(c)) = ...
            out.TestGLMs{subjChans(c)}.Rsquared.Adjusted;
    end
end

%% K-Means on the + and - components
out.GLMsCompsKMeans     = 3;
out.KMeansReplicates    = 100;

%% Study
out.StudyGLMsCompKmeans =[];
X = out.StudyGLMSChanCompTVal;

% Study Positive
[ch,co] = find(X>0);
out.StudyGLMsCompKmeans.PosCompIDs = [ch,co];
nPosComps = numel(ch);
Y = zeros(nPosComps,out.nFeat);
for ii = 1:nPosComps
    Y(ii,:) = out.Projections(ch(ii),:,co(ii));
end
[IDX, C, SUMD, D] =...
    kmeans(Y, out.GLMsCompsKMeans,'replicates',out.KMeansReplicates,'distance','correlation');
out.StudyGLMsCompKmeans.PIDX     = IDX;
out.StudyGLMsCompKmeans.PC       = C;
out.StudyGLMsCompKmeans.PSUMD    = SUMD;
out.StudyGLMsCompKmeans.PD       = D;

% Get Tvals for every cluster
for ii = 1:out.GLMsCompsKMeans
    clcomps = find(IDX==ii); % cluster IDs for each componenent
    clch = ch(clcomps) ;
    clco= co(clcomps);
    nn = numel(clcomps);
    % create matrix for each cluster
    Z = zeros(nn,out.nFeat);
    for jj = 1:nn
        Z(jj,:) = out.Projections(clch(jj),:,clco(jj));
    end
    [~,p,~,t] = ttest(Z);
    out.StudyGLMsCompKmeans.PT(ii,:) = t.tstat;
    out.StudyGLMsCompKmeans.PP(ii,:) = p;
end

% Study Negative
[ch,co] = find(X<0);
out.StudyGLMsCompKmeans.NegCompIDs = [ch,co];
nNegComps = numel(ch);
Y = zeros(nNegComps,out.nFeat);
for ii = 1:nNegComps
    Y(ii,:) = out.Projections(ch(ii),:,co(ii));
end
[IDX, C, SUMD, D] =...
    kmeans(Y, out.GLMsCompsKMeans,'replicates',out.KMeansReplicates,'distance','correlation');
out.StudyGLMsCompKmeans.NIDX     = IDX;
out.StudyGLMsCompKmeans.NC       = C;
out.StudyGLMsCompKmeans.NSUMD    = SUMD;
out.StudyGLMsCompKmeans.ND       = D;

% Get Tvals for every cluster
for ii = 1:out.GLMsCompsKMeans
    clcomps = find(IDX==ii); % cluster IDs for each componenent
    clch = ch(clcomps) ;
    clco= co(clcomps);
    nn = numel(clcomps);
    % create matrix for each cluster
    Z = zeros(nn,out.nFeat);
    for jj = 1:nn
        Z(jj,:) = out.Projections(clch(jj),:,clco(jj));
    end
    [~,p,~,t] = ttest(Z);
    out.StudyGLMsCompKmeans.NT(ii,:) = t.tstat;
    out.StudyGLMsCompKmeans.NP(ii,:) = p;
end

%% Test
out.TestGLMsCompKmeans =[];
X = out.TestGLMSChanCompTVal;
% Test Positive
[ch,co] = find(X>out.GLMsCompsThr);
out.TestGLMsCompKmeans.PosCompIDs = [ch,co];
nPosComps = numel(ch);
Y = zeros(nPosComps,out.nFeat);
for ii = 1:nPosComps
    Y(ii,:) = out.Projections(ch(ii),:,co(ii));
end
[IDX, C, SUMD, D] =...
    kmeans(Y, out.GLMsCompsKMeans,'replicates',out.KMeansReplicates,'distance','correlation');
out.TestGLMsCompKmeans.PIDX     = IDX;
out.TestGLMsCompKmeans.PC       = C;
out.TestGLMsCompKmeans.PSUMD    = SUMD;
out.TestGLMsCompKmeans.PD       = D;
% Get Tvals for every cluster
for ii = 1:out.GLMsCompsKMeans
    clcomps = find(IDX==ii); % cluster IDs for each componenent
    clch = ch(clcomps) ;
    clco= co(clcomps);
    nn = numel(clcomps);
    % create matrix for each cluster
    Z = zeros(nn,out.nFeat);
    for jj = 1:nn
        Z(jj,:) = out.Projections(clch(jj),:,clco(jj));
    end
    [~,p,~,t] = ttest(Z);
    out.TestGLMsCompKmeans.PT(ii,:) = t.tstat;
    out.TestGLMsCompKmeans.PP(ii,:) = p;
end

% Test Negative
[ch,co] = find(X<-out.GLMsCompsThr);
out.TestGLMsCompKmeans.NegCompIDs = [ch,co];
nNegComps = numel(ch);
Y = zeros(nNegComps,out.nFeat);
for ii = 1:nNegComps
    Y(ii,:) = out.Projections(ch(ii),:,co(ii));
end
[IDX, C, SUMD, D] =...
    kmeans(Y, out.GLMsCompsKMeans,'replicates',out.KMeansReplicates,'distance','correlation');
out.TestGLMsCompKmeans.NIDX      = IDX;
out.TestGLMsCompKmeans.NC       = C;
out.TestGLMsCompKmeans.NSUMD    = SUMD;
out.TestGLMsCompKmeans.ND       = D;
% Get Tvals for every cluster
for ii = 1:out.GLMsCompsKMeans
    clcomps = find(IDX==ii); % cluster IDs for each componenent
    clch = ch(clcomps) ;
    clco= co(clcomps);
    nn = numel(clcomps);
    % create matrix for each cluster
    Z = zeros(nn,out.nFeat);
    for jj = 1:nn
        Z(jj,:) = out.Projections(clch(jj),:,clco(jj));
    end
    [~,p,~,t] = ttest(Z);
    out.TestGLMsCompKmeans.NT(ii,:) = t.tstat;
    out.TestGLMsCompKmeans.NP(ii,:) = p;
end

%% Get Distribution of Selected components by ROIs
K       = out.GLMsCompsKMeans;
rois    = data.ROIid;
nROIs   = numel(unique(rois));

StudTest    = {'StudyGLMsCompKmeans','TestGLMsCompKmeans'};
PosNegComp  = {'PosCompIDs','NegCompIDs'};
PosNegCLID  = {'PIDX','NIDX'};

% components above/below threshold
comps = cell(2,2); % rows: pos/neg: col: study/test
out.GLMsThrComps  = [];
for ii =1:2
    for jj=1:2
        comps{ii,jj}=out.(StudTest{ii}).(PosNegComp{jj})(:,1);
    end
end
out.GLMsThrComps.CompsStr = {'+Study','+Test';'-Study','-Test'};
out.GLMsThrComps.CompIDs = comps;

% Distribution of rois by pos/neg-study/test components.
Ns = zeros(nROIs,2,2);
% Distribution of rois by cluster AND pos/neg-study/test components.
ROIkCompNs = zeros(nROIs,K,2,2);
for rr=1:nROIs
    for ii = 1:2
        for jj = 1:2
            for kk=1:K
                cLIDs = out.(StudTest{ii}).(PosNegCLID{jj})==kk;
                ROIkCompNs(rr,kk,ii,jj) = sum(rois(comps{ii,jj}(cLIDs))==rr);
            end % clusters
            Ns(rr,ii,jj)=sum(rois(comps{ii,jj})==rr);
        end % study/test
    end % pos/neg
end % rois
out.GLMsThrComps.nCompsByCLROI = ROIkCompNs;
out.GLMsThrComps.nCompsByCLROIDims = {'ROI','K','Study/Test','Pos/Neg'};


ROIKCont = zeros(nROIs,K,2,2);
for kk = 1:K
    temp    = squeeze(ROIkCompNs(:,kk,:,:))./Ns;
    tempS   = squeeze(sum(temp));
    for ii =1:2
        for jj=1:2
            ROIKCont(:,kk,ii,jj) = temp(:,ii,jj)/tempS(ii,jj);
        end
    end
end
out.GLMsThrComps.relContCLROI = ROIKCont;



