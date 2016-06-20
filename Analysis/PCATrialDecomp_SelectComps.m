function out = PCATrialDecomp_SelectComps(data,opts)
% this analysis function finds the principal components for each channel
% that explains the most variance across trials.
%
% it also uses these components to
%
% opts:
% 	nComps -> number of components
%
% data input must be from the groupLPCDataMultiband.

out                 = [];
out.nSubjs 			= size(data.BinERPs,1);
out.nComps 			= opts.nComps;
out.nChans 			= data.nChans;
out.nFeat           = data.nBands*sum(data.AnalysisBins);
out.ROIs            = data.ROIid; rois =data.ROIid;
out.nROIs           = numel(unique(data.ROIid));

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
out.StudyGLMsChanR2         = nan(out.nChans,1);
out.StudyGLMsChanAR2         = nan(out.nChans,1);

out.TestSelectedComps       = cell(out.nChans,1);
out.TestGLMs        		= cell(out.nChans,1);
out.TestGLMSChanCompTVal 	= nan(out.nChans,out.nComps);
out.TestGLMsChanR2          = nan(out.nChans,1);
out.TestGLMsChanAR2         = nan(out.nChans,1);

for ss=1:out.nSubjs
    subjChans = find(data.subjChans==ss);
    % re-order to channels, trials , bands , time bins
    x = permute(data.BinERPs{ss}(:,:,:,data.AnalysisBins),[2 3 1 4]);
    % concatenate bands and time bins
    x = x(:,:,:);
    [nSubjChan,nSubjTrials,~] = size(x); % d1 -> n subj channels, d2 -> # of trials, d3 -> bands*timebins
    out.sRTs{ss} = data.studyRTs{ss}(data.trials{ss});
    out.tRTs{ss} = data.testRTs{ss}(data.trials{ss});
    rts1 = -log10(out.sRTs{ss});
    rts2 = -log10(out.tRTs{ss});
    
    out.Comps{ss} = zeros(nSubjChan,nSubjTrials,out.nComps);
    
    
    % for subject channels
    for c=1:nSubjChan
        % get the trial by band*bin matrix for each channel.
        A = squeeze(x(c,:,:));
        
        % PCA
        [C,S,~,~,E]= pca(A','NumComponents',out.nComps,'Centered',false);
        % correct for sign (based on inner product of residuals
        [C,S] = getCompSign(C,S);
        
        out.Comps{ss}(c,:,:) = C;
        out.Projections(subjChans(c),:,:) = S;
        out.VarExp(subjChans(c),:) = E(1:out.nComps);
        
        % Correlations
        [r1,p1] = corr(C,rts1,'type','spearman');
        [r2,p2] = corr(C,rts2,'type','spearman');
        out.CorrStudyRTs(subjChans(c),:) = r1;
        out.CorrTestRTs(subjChans(c),:)  = r2;
        out.CorrStudyRTsP(subjChans(c),:) = p1;
        out.CorrTestRTsP(subjChans(c),:)  = p2;
        
        CompsIdx1 = find(abs(r1)>=opts.rThr);
        out.StudySelectedComps{subjChans(c)} = CompsIdx1;
        CompsIdx2 = find(abs(r2)>=opts.rThr);
        out.TestSelectedComps{subjChans(c)} = CompsIdx2;
        
        % GLMSs
        if ~isempty(CompsIdx1)
            out.StudyGLMs{subjChans(c)} = fitglm(C(:,CompsIdx1),rts1);
            out.StudyGLMSChanCompTVal(subjChans(c),CompsIdx1) ...
                = out.StudyGLMs{subjChans(c)}.Coefficients.tStat(2:end);
            out.StudyGLMsChanR2(subjChans(c))...
                = out.StudyGLMs{subjChans(c)}.Rsquared.Ordinary;
            out.StudyGLMsChanAR2(subjChans(c)) = ...
                out.StudyGLMs{subjChans(c)}.Rsquared.Adjusted;
        end
        if ~isempty(CompsIdx2)
            out.TestGLMs{subjChans(c)} = fitglm(C(:,CompsIdx2),rts2);
            out.TestGLMSChanCompTVal(subjChans(c),CompsIdx2) ...
                = out.TestGLMs{subjChans(c)}.Coefficients.tStat(2:end);
            out.TestGLMsChanR2(subjChans(c)) ...
                = out.TestGLMs{subjChans(c)}.Rsquared.Ordinary;
            out.TestGLMsChanAR2(subjChans(c)) = ...
                out.TestGLMs{subjChans(c)}.Rsquared.Adjusted;
        end
    end
end

%% Matrix of Components that explain RTs: Study
out.StudySelComps =[];
X = out.CorrStudyRTs;

[chP,coP] = find(X>opts.rThr);
out.StudySelComps.PosCompIDs = [chP,coP];
[chN,coN] = find(X<-opts.rThr );
out.StudySelComps.NegCompIDs = [chN,coN];
ch = [chP;chN]; co = [coP;coN];
out.StudySelComps.CompIDs    = [ch,co];
nPosComps = numel(chP);
nNegComps = numel(chN);
nComps  = nPosComps+nNegComps;

Y = zeros(nComps,out.nFeat);
out.StudySelComps.Corrs = zeros(nComps,1);
for ii = 1:nComps
    if ii <= nPosComps
        Y(ii,:) = out.Projections(ch(ii),:,co(ii));
        out.StudySelComps.Corrs(ii) = X(ch(ii),co(ii));
    else
        Y(ii,:) = -out.Projections(ch(ii),:,co(ii));
        out.StudySelComps.Corrs(ii) = -X(ch(ii),co(ii));
    end
    out.StudySelComps.Corrs(ii) = X(ch(ii),co(ii));
end
Yp = Y(1:nPosComps,:);
Xp = out.StudySelComps.Corrs(1:nPosComps,:);
[c,p]=corr(Yp,Xp);
out.StudySelComps.PCompCorr = [c,p];

Yn = Y((nPosComps+1):end,:);
Xn = out.StudySelComps.Corrs((nPosComps+1):end,:);
[c,p]=corr(Yn,Xn);
out.StudySelComps.NCompCorr = [c,p];

Y2 = [Yp;-Yn]; X2 = [Xp;-Xn];
[c,p]=corr(Y2,X2);
out.StudySelComps.CompCorr = [c,p];

out.StudySelComps.Mat = Y;
%% %% Find Covariance Among Selected Components by Region.

%components across all selected components
out.StudyPCASelComps =[];
[C,S,~,~,E]= pca(Y','NumComponents',out.nComps,'centered',false);
[C,S] = getCompSign(C,S);
out.StudyPCASelComps.C =C;
out.StudyPCASelComps.S =S;
out.StudyPCASelComps.E =E;

% components by region
out.StudyPCASelCompsROIs =[];
out.StudyPCASelCompsROIs.C =cell(out.nROIs,1);
out.StudyPCASelCompsROIs.S =cell(out.nROIs,1);
out.StudyPCASelCompsROIs.E =cell(out.nROIs,1);
out.StudyPCASelCompsROIs.compROIs =cell(out.nROIs,1);
for rr = 1:out.nROIs
    compROIs = find(rois(ch)==rr);
    [C,S,~,~,E]= pca(Y(compROIs,:)','NumComponents',out.nComps,'centered',false);
    [C,S] = getCompSign(C,S);
    out.StudyPCASelCompsROIs.C{rr} = C;
    out.StudyPCASelCompsROIs.S{rr} = S;
    out.StudyPCASelCompsROIs.E{rr} = E;
    out.StudyPCASelCompsROIs.compROIs{rr} = compROIs;
        
    pCompROIs = compROIs(compROIs<=nPosComps);
    Yp = Y(pCompROIs,:);
    Xp = out.StudySelComps.Corrs(pCompROIs,:);
    [c,p]=corr(Yp,Xp);
    out.StudyPCASelCompsROIs.PCompCorr{rr} = [c,p];
    
    nCompROIs = compROIs(compROIs>nPosComps);
    Yn = Y(nCompROIs,:);
    Xn = out.StudySelComps.Corrs(nCompROIs,:);
    [c,p]=corr(Yn,Xn);
    out.StudyPCASelCompsROIs.NCompCorr{rr} = [c,p];
    
    Y2 = [Yp;-Yn]; X2 = [Xp;-Xn];
    [c,p]=corr(Y2,X2);
    out.StudyPCASelCompsROIs.CompCorr{rr} = [c,p];
end


%% Matrix of Components that explain RTs: Test
out.TestSelComps =[];
X = out.CorrTestRTs;

% Test Positive
[chP,coP] = find(X>opts.rThr);
out.TestSelComps.PosCompIDs = [chP,coP];
[chN,coN] = find(X<-opts.rThr);
out.TestSelComps.NegCompIDs = [chN,coN];
ch = [chP;chN]; co = [coP;coN];
out.TestSelComps.CompIDs    = [ch,co];

nPosComps = numel(chP);
nNegComps = numel(chN);
nComps  = nPosComps+nNegComps;

Y  = zeros(nComps,out.nFeat);
out.TestSelComps.Corrs = zeros(nComps,1);
for ii = 1:nComps
    if ii <= nPosComps
        Y(ii,:) = out.Projections(ch(ii),:,co(ii));
        out.TestSelComps.Corrs(ii) = X(ch(ii),co(ii));
    else
        Y(ii,:) = -out.Projections(ch(ii),:,co(ii));
        out.TestSelComps.Corrs(ii) = -X(ch(ii),co(ii));
    end
    
end
Yp = Y(1:nPosComps,:);
Xp = out.TestSelComps.Corrs(1:nPosComps,:);
[c,p]=corr(Yp,Xp);
out.TestSelComps.PCompCorr = [c,p];

Yn = Y((nPosComps+1):end,:);
Xn = out.TestSelComps.Corrs((nPosComps+1):end,:);
[c,p]=corr(Yn,Xn);
out.TestSelComps.NCompCorr = [c,p];

Y2 = [Yp;-Yn]; X2 = [Xp;-Xn];
[c,p]=corr(Y2,X2);
out.TestSelComps.CompCorr = [c,p];
out.TestSelComps.Mat = Y;

%% Find Covariance Among Selected Components by Region.
%components across all selected components
out.TestPCASelComps =[];
[C,S,~,~,E]= pca(Y','NumComponents',out.nComps,'centered',false);
out.TestPCASelComps.Corig =C;
out.TestPCASelComps.Sorig =S;
[C,S] = getCompSign(C,S);
out.TestPCASelComps.C =C;
out.TestPCASelComps.S =S;
out.TestPCASelComps.E =E;

% components by region
out.TestPCASelCompsROIs =[];
out.TestPCASelCompsROIs.C =cell(out.nROIs,1);
out.TestPCASelCompsROIs.S =cell(out.nROIs,1);
out.TestPCASelCompsROIs.E =cell(out.nROIs,1);
out.TestPCASelCompsROIs.compROIs = cell(out.nROIs,1);
for rr = 1:out.nROIs
    compROIs = find(rois(ch)==rr);
    [C,S,~,~,E]= pca(Y(compROIs,:)','NumComponents',out.nComps,'centered',false);
    [C,S] = getCompSign(C,S);
    out.TestPCASelCompsROIs.C{rr} = C;
    out.TestPCASelCompsROIs.S{rr} = S;
    out.TestPCASelCompsROIs.E{rr} = E;
    out.TestPCASelCompsROIs.compROIs{rr} = compROIs;
    
    pCompROIs = compROIs(compROIs<=nPosComps);
    Yp = Y(pCompROIs,:);
    Xp = out.TestSelComps.Corrs(pCompROIs,:);
    [c,p]=corr(Yp,Xp);
    out.TestPCASelCompsROIs.PCompCorr{rr} = [c,p];
    
    nCompROIs = compROIs(compROIs>nPosComps);
    Yn = Y(nCompROIs,:);
    Xn = out.TestSelComps.Corrs(nCompROIs,:);
    [c,p]=corr(Yn,Xn);    
    out.TestPCASelCompsROIs.NCompCorr{rr} = [c,p];    
    
    Y2 = [Yp;-Yn]; X2 = [Xp;-Xn];
    [c,p]=corr(Y2,X2);
    out.TestPCASelCompsROIs.CompCorr{rr} = [c,p];
end

end
function [C,S]=getCompSign(C,S)

nComps = size(C,2);
try
    for cc = 1:nComps
        
        x=mean(C(:,cc)*S(:,cc)')';
        negComp = sign(corr(x,S(:,cc)));
        C(:,cc) = negComp*C(:,cc);
        S(:,cc) = x;
        %     Y=X-C(:,setdiff(1:nComps,cc))*S(:,setdiff(1:nComps,cc))';
        %     compSgn=sign(mean(Y*S(:,cc)));
        %     S(:,cc) = compSgn*S(:,cc);
        %     C(:,cc) = compSgn*C(:,cc);
        %
        %     % Get correct scaling for componenets
        %     S(:,cc) = mean(X*S(:,cc)')';
    end
 
    
catch
end
end
