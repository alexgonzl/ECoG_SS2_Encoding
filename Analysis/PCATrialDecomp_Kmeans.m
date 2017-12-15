function out = PCATrialDecomp_Kmeans(out,opts)

% set random seed for stability of clusters.
rng(1);
K                       = opts.maxK;
out.GLMsCompsKMeans     = K  ;
out.KMeansReplicates    = 100;
out.KSelThr             = 0.05; % increment ncecessary

rois = out.ROIs; nROIs = out.nROIs;
%% study: Kmeans for all selected components

Y = out.StudySelComps.Mat;

% Get OptimalK
[KS,out.StudyKmeansScores , out.StudyKmeansScoresIncr] = getOptimalK(Y,K,out.KSelThr);
out.StudyKmeanOptK = KS;

fprintf('Optimal K for Study: %i \n',KS)
out.StudyCompKmeans = [];
[IDX, C, SUMD, D ]=...
    kmeans(Y, KS,'replicates',out.KMeansReplicates,'distance','correlation');
out.StudyCompKmeans.IDX     = IDX;
out.StudyCompKmeans.C       = C;
out.StudyCompKmeans.SUMD    = SUMD;
out.StudyCompKmeans.D       = D;
out.StudyCompKmeans.CompClCorrs = cell(KS,1);
for ii = 1:KS
    clcomps = find(IDX==ii);
    out.StudyCompKmeans.CompClCorrs{ii} = out.StudySelComps.Corrs(clcomps);
    Z = Y(clcomps,:);
    [~,p,~,t] = ttest(Z);
    out.StudyCompKmeans.T(ii,:) = t.tstat;
    out.StudyCompKmeans.P(ii,:) = p;
end

%% study: Kmeans by ROI

out.StudyKmeanCompROIs      = [];
out.StudyKmeanCompROIs.Ko   = nan(nROIs,1);
out.StudyKmeanCompROIs.CL_Ns= cell(nROIs,1);
out.StudyKmeanCompROIs.IDX  = cell(nROIs,1);
out.StudyKmeanCompROIs.C    = cell(nROIs,1);
out.StudyKmeanCompROIs.SUMD = cell(nROIs,1);
out.StudyKmeanCompROIs.D    = cell(nROIs,1);
out.StudyKmeanCompROIs.compROIs = cell(nROIs,1);
out.StudyKmeanCompROIs.T    = cell(nROIs,1);
out.StudyKmeanCompROIs.P    = cell(nROIs,1);
out.StudyKmeanCompROIs.CompClCorrs = [];
chanIDs = out.StudySelComps.CompIDs(:,1);

for rr = 1:nROIs
    compROIs = find(rois(chanIDs)==rr);
    out.StudyKmeanCompROIs.compROIs{rr} = compROIs;
    Yr = Y(compROIs,:);
    
    % Get OptimalK per region
    Kr = getOptimalK(Yr,K,out.KSelThr);
    fprintf('Optimal K for ROI %i Study: %i \n',rr,Kr)
    out.StudyKmeanCompROIs.Ko(rr) = Kr;
    [IDX, C, SUMD, D ]=...
        kmeans(Yr, Kr,'replicates',out.KMeansReplicates,'distance','correlation');
    
    out.StudyKmeanCompROIs.IDX{rr}     = IDX;
    out.StudyKmeanCompROIs.CL_Ns{rr}   = histc(IDX,1:Kr);
    out.StudyKmeanCompROIs.C{rr}       = C;
    out.StudyKmeanCompROIs.SUMD{rr}    = SUMD;
    out.StudyKmeanCompROIs.D{rr}       = D;
    for kk = 1:Kr
        roiclcomps = compROIs(IDX==kk);
        out.StudyKmeanCompROIs.CompClCorrs{rr,kk} = out.StudySelComps.Corrs(roiclcomps);
        Z = Y(roiclcomps,:);
        [~,p,~,t] = ttest(Z);
        out.StudyKmeanCompROIs.T{rr}(kk,:) = t.tstat;
        out.StudyKmeanCompROIs.P{rr}(kk,:) = p;
    end
end
%% test: Kmeans for all selected components

Y = out.TestSelComps.Mat;

% Get OptimalK
[KT,out.TestKmeansScores , out.TestKmeansScoresIncr] = getOptimalK(Y,K,out.KSelThr);
out.TestKmeanOptK = KT;

fprintf('Optimal K for Test: %i \n',KT)
out.TestCompKmeans = [];
[IDX, C, SUMD, D ]=...
    kmeans(Y, KT,'replicates',out.KMeansReplicates,'distance','correlation');

out.TestCompKmeans.IDX     = IDX;
out.TestCompKmeans.C       = C;
out.TestCompKmeans.SUMD    = SUMD;
out.TestCompKmeans.D       = D;
out.TestCompKmeans.CompClCorrs = cell(KT,1);

for ii = 1:KT
    clcomps = find(IDX==ii);
    out.TestCompKmeans.CompClCorrs{ii} = out.TestSelComps.Corrs(clcomps);
    Z = Y(clcomps,:);
    [~,p,~,t] = ttest(Z);
    out.TestCompKmeans.T(ii,:) = t.tstat;
    out.TestCompKmeans.P(ii,:) = p;
end

%% Test: Kmeans by ROI

out.TestKmeanCompROIs      = [];
out.TestKmeanCompROIs.Ko   = nan(nROIs,1);
out.TestKmeanCompROIs.CL_Ns= cell(nROIs,1);
out.TestKmeanCompROIs.IDX  = cell(nROIs,1);
out.TestKmeanCompROIs.C    = cell(nROIs,1);
out.TestKmeanCompROIs.SUMD = cell(nROIs,1);
out.TestKmeanCompROIs.D    = cell(nROIs,1);
out.TestKmeanCompROIs.compROIs = cell(nROIs,1);
out.TestKmeanCompROIs.T    = cell(nROIs,1);
out.TestKmeanCompROIs.P    = cell(nROIs,1);

chanIDs = out.TestSelComps.CompIDs(:,1);

for rr = 1:nROIs
    compROIs = find(rois(chanIDs)==rr);
    out.TestKmeanCompROIs.compROIs{rr} = compROIs;
    Yr = Y(compROIs,:);
    
    % Get OptimalK per region
    
    Kr = getOptimalK(Yr,K,out.KSelThr);
    if isnan(Kr)
        Kr=5;
    end
    out.TestKmeanCompROIs.Ko(rr) = Kr;
    fprintf('Optimal K for ROI %i Tesrt: %i \n',rr,Kr)
    [IDX, C, SUMD, D ]=...
        kmeans(Yr, Kr,'replicates',out.KMeansReplicates,'distance','correlation');
    
    out.TestKmeanCompROIs.IDX{rr}     = IDX;
    out.TestKmeanCompROIs.CL_Ns{rr}   = histc(IDX,1:Kr);
    out.TestKmeanCompROIs.C{rr}       = C;
    out.TestKmeanCompROIs.SUMD{rr}    = SUMD;
    out.TestKmeanCompROIs.D{rr}       = D;
    for kk = 1:Kr
        roiclcomps = compROIs(IDX==kk);
        out.TestKmeanCompROIs.CompClCorrs{rr,kk} = out.TestSelComps.Corrs(roiclcomps);
        Z = Y(roiclcomps,:);
        [~,p,~,t] = ttest(Z);
        out.TestKmeanCompROIs.T{rr}(kk,:) = t.tstat;
        out.TestKmeanCompROIs.P{rr}(kk,:) = p;
    end
end
end

function [Ko,Scores,ScoresInc] = getOptimalK(Y,maxK,Kthr)
Ko          = nan;
Scores      = nan(maxK,1);
ScoresInc   = nan(maxK,1);

for kk = 1:maxK
    [IDX, ~, SUMD, ~] =...
        kmeans(Y, kk,'replicates',100,'distance','correlation');
    Ns = histc(IDX,1:kk);
    Scores(kk) = mean(SUMD./Ns);
    if kk>2
        ScoresInc(kk) = (Scores(kk-1)-Scores(kk))/Scores(kk-1);
        if ScoresInc(kk)<Kthr
            Ko= kk-1;
            break
        end
    end
end
return
end
