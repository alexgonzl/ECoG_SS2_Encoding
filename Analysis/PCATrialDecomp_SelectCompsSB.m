function out = PCATrialDecomp_SelectCompsSB(data,opts)
% this analysis function finds the principal components for each channel
% that explains the most variance across trials, for each band
%
% data input must be from the groupLPCDataMultiband.

out                 = [];
out.nSubjs 			= size(data.BinERPs,1);
out.nComps 			= opts.nComps;
out.rThr            = opts.rThr;
out.nChans 			= data.nChans;
out.nBands          = data.nBands;
out.nFeat           = sum(data.AnalysisBins);
out.ROIs            = data.ROIid; rois =data.ROIid;
out.nROIs           = numel(unique(data.ROIid));

% Pre allocation for outputs from PCA
out.Comps           = cell(out.nSubjs);
out.Projections     = zeros(out.nChans,out.nBands,out.nFeat,out.nComps);
out.VarExp          = zeros(out.nChans,out.nBands,out.nComps);

% Pre allocation for correlations
out.CorrStudyRTs    = zeros(out.nChans,out.nBands,out.nComps);
out.CorrTestRTs     = zeros(out.nChans,out.nBands,out.nComps);

out.StudySelectedComps       = cell(out.nChans,out.nBands);
out.StudyGLMs        		= cell(out.nChans,out.nBands);
out.StudyGLMSChanCompTVal 	= nan(out.nChans,out.nBands,out.nComps);
out.StudyGLMsChanR2          = nan(out.nChans,out.nBands);
out.StudyGLMsChanAR2         = nan(out.nChans,out.nBands);

out.TestSelectedComps       = cell(out.nChans,out.nBands);
out.TestGLMs        		= cell(out.nChans,out.nBands);
out.TestGLMSChanCompTVal 	= nan(out.nChans,out.nBands,out.nComps);
out.TestGLMsChanR2          = nan(out.nChans,out.nBands);
out.TestGLMsChanAR2         = nan(out.nChans,out.nBands);


for ss=1:out.nSubjs
    subjChans = find(data.subjChans==ss);
    % re-order to channels, trials , bands , time bins
    x = permute(data.BinERPs{ss}(:,:,:,data.AnalysisBins),[2 3 1 4]);
    
    [nSubjChan,nSubjTrials,~,~] = size(x); % d1 -> n subj channels, d2 -> # of trials, d3 -> bands*timebins
    rts1 = -log10(data.studyRTs{ss}(data.trials{ss}));
    rts2 = -log10(data.testRTs{ss}(data.trials{ss}));
    
    out.Comps{ss} = zeros(nSubjChan,out.nBands,nSubjTrials,out.nComps);
    
    
    for bb = 1:out.nBands
    % for subject channels
    for c=1:nSubjChan
        % get the trial by band*bin matrix for each channel.
        A = squeeze(x(c,:,bb,:));
        
        % PCA
        [C,S,~,~,E]= pca(A','NumComponents',out.nComps,'Centered',false);
        % correct for sign (based on inner product of residuals
        [C,S] = getCompSign(C,S);
        
        out.Comps{ss}(c,bb,:,:) = C;
        out.Projections(subjChans(c),bb,:,:) = S;
        out.VarExp(subjChans(c),bb,:) = E(1:out.nComps);
        
        % Correlations
        r1 = corr(C,rts1,'type','spearman');
        r2 = corr(C,rts2,'type','spearman');
        out.CorrStudyRTs(subjChans(c),bb,:) = r1;
        out.CorrTestRTs(subjChans(c),bb,:)  = r2;
        
        CompsIdx1 = find(abs(r1)>=out.rThr);
        out.StudySelectedComps{subjChans(c),bb} = CompsIdx1;
        CompsIdx2 = find(abs(r2)>=out.rThr);
        out.TestSelectedComps{subjChans(c),bb} = CompsIdx2;
        
        % GLMSs
        if ~isempty(CompsIdx1)
            out.StudyGLMs{subjChans(c),bb} = fitglm(C(:,CompsIdx1),rts1);
            out.StudyGLMSChanCompTVal(subjChans(c),bb,CompsIdx1) ...
                = out.StudyGLMs{subjChans(c),bb}.Coefficients.tStat(2:end);
            out.StudyGLMsChanR2(subjChans(c),bb)...
                = out.StudyGLMs{subjChans(c),bb}.Rsquared.Ordinary;
            out.StudyGLMsChanAR2(subjChans(c),bb) = ...
                out.StudyGLMs{subjChans(c),bb}.Rsquared.Adjusted;
        end
        if ~isempty(CompsIdx2)
            out.TestGLMs{subjChans(c),bb} = fitglm(C(:,CompsIdx2),rts2);
            out.TestGLMSChanCompTVal(subjChans(c),bb,CompsIdx2) ...
                = out.TestGLMs{subjChans(c),bb}.Coefficients.tStat(2:end);
            out.TestGLMsChanR2(subjChans(c),bb) ...
                = out.TestGLMs{subjChans(c),bb}.Rsquared.Ordinary;
            out.TestGLMsChanAR2(subjChans(c),bb) = ...
                out.TestGLMs{subjChans(c),bb}.Rsquared.Adjusted;
        end
    end
    end
end

%% Matrix of Components that explain RTs: Study
out.StudySelComps =[];
for bb = 1:out.nBands
X = squeeze(out.CorrStudyRTs(:,bb,:));

[chP,coP] = find(X>out.rThr);
out.StudySelComps.PosCompIDs{bb} = [chP,coP];
[chN,coN] = find(X<-out.rThr );
out.StudySelComps.NegCompIDs{bb} = [chN,coN];
ch = [chP;chN]; co = [coP;coN];
out.StudySelComps.CompIDs{bb}    = [ch,co];
nPosComps = numel(chP);
nNegComps = numel(chN);
nComps  = nPosComps+nNegComps;

Y = zeros(nComps,out.nFeat);
out.StudySelComps.Corrs{bb} = nan(nComps,1);
for ii = 1:nComps
    if ii <= nPosComps
        Y(ii,:) = out.Projections(ch(ii),bb,:,co(ii));
    else
        Y(ii,:) = -out.Projections(ch(ii),bb,:,co(ii));        
    end
    out.StudySelComps.Corrs{bb}(ii) = X(ch(ii),co(ii));
end
Yp = Y(1:nPosComps,:);
Xp = out.StudySelComps.Corrs{bb}(1:nPosComps,:);
[c,p]=corr(Yp,Xp);
out.StudySelComps.PCompCorr{bb} = [c,p];

Yn = Y((nPosComps+1):end,:);
Xn = out.StudySelComps.Corrs{bb}((nPosComps+1):end,:);
[c,p]=corr(Yn,Xn);
out.StudySelComps.NCompCorr{bb} = [c,p];

Y2 = [Yp;-Yn]; X2 = [Xp;Xn];
[c,p]=corr(Y2,X2);
out.StudySelComps.CompCorr{bb} = [c,p];

out.StudySelComps.Mat{bb} = Y;
end
for bb=1:out.nBands    
    [~,p,~,t]=ttest(out.StudySelComps.Mat{bb});
    out.StudySelComps.T(bb,:) = t.tstat;
    out.StudySelComps.P(bb,:) = p;
    ch = out.StudySelComps.CompIDs{bb}(:,1);
    for rr = 1:out.nROIs        
        compROIs = find(rois(ch)==rr);
        Z=out.StudySelComps.Mat{bb}(compROIs,:);
        [~,p,~,t]=ttest(Z);
        out.StudySelComps.TR(rr,bb,:) = t.tstat;
        out.StudySelComps.PR(rr,bb,:) = p;        
        [C,S,~,~,E]= pca(Z','NumComponents',out.nComps,'Centered',false);
        [C,S] = getCompSign(C,S);
        out.StudySelComps.PCs_S(rr,bb,:,:) = (S);
        out.StudySelComps.PCs_C{rr,bb} = C;
        out.StudySelComps.PCs_E(rr,bb,:)   = E(1:out.nComps);
    end
end

%% Matrix of Components that explain RTs: Test
out.TestSelComps =[];

for bb=1:out.nBands
X = squeeze(out.CorrTestRTs(:,bb,:));

% Test Positive
[chP,coP] = find(X>opts.rThr);
out.TestSelComps.PosCompIDs{bb} = [chP,coP];
[chN,coN] = find(X<-opts.rThr);
out.TestSelComps.NegCompIDs{bb} = [chN,coN];
ch = [chP;chN]; co = [coP;coN];
out.TestSelComps.CompIDs{bb}    = [ch,co];

nPosComps = numel(chP);
nNegComps = numel(chN);
nComps  = nPosComps+nNegComps;

Y  = zeros(nComps,out.nFeat);
out.TestSelComps.Corrs{bb} = zeros(nComps,1);
for ii = 1:nComps
    if ii <= nPosComps
        Y(ii,:) = out.Projections(ch(ii),bb,:,co(ii));        
    else
        Y(ii,:) = -out.Projections(ch(ii),bb,:,co(ii));
    end
    out.TestSelComps.Corrs{bb}(ii) = X(ch(ii),co(ii));
end
Yp = Y(1:nPosComps,:);
Xp = out.TestSelComps.Corrs{bb}(1:nPosComps,:);
[c,p]=corr(Yp,Xp);
out.TestSelComps.PCompCorr{bb} = [c,p];

Yn = Y((nPosComps+1):end,:);
Xn = out.TestSelComps.Corrs{bb}((nPosComps+1):end,:);
[c,p]=corr(Yn,-Xn);
out.TestSelComps.NCompCorr{bb} = [c,p];

Y2 = [Yp;-Yn]; X2 = [Xp;Xn];
[c,p]=corr(Y2,X2);
out.TestSelComps.CompCorr = [c,p];
out.TestSelComps.Mat{bb} = Y;
end
for bb=1:out.nBands
    [~,p,~,t]=ttest(out.TestSelComps.Mat{bb});
    out.TestSelComps.T(bb,:) = t.tstat;
    out.TestSelComps.P(bb,:) = p;
    ch = out.TestSelComps.CompIDs{bb}(:,1);
    for rr = 1:out.nROIs        
        compROIs = find(rois(ch)==rr);
         Z=out.TestSelComps.Mat{bb}(compROIs,:);
        [~,p,~,t]=ttest(Z);        
        out.TestSelComps.TR(rr,bb,:) = t.tstat;
        out.TestSelComps.PR(rr,bb,:) = p;
        
        [C,S,~,~,E]= pca(Z','NumComponents',out.nComps,'Centered',false);
        [C,S] = getCompSign(C,S);
        out.TestSelComps.PCs_S(rr,bb,:,:) = (S);
        out.TestSelComps.PCs_C{rr,bb} = C;
        out.TestSelComps.PCs_E(rr,bb,:)   = E(1:out.nComps);
    end    
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
