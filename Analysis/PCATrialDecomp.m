function out = PCATrialDecomp(data,opts)
% this analysis function finds the principal components for each channel
% that explains the most variance across channels.
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

% Pre allocation for outputs from PCA
out.Comps           = cell(out.nSubjs,1);
out.Projections     = zeros(out.nChans,out.nFeat,out.nComps);
out.VarExp          = zeros(out.nChans,out.nComps);

% Pre allocation for correlations
out.CorrStudyRTs    = zeros(out.nChans,out.nComps);
out.CorrTestRTs     = zeros(out.nChans,out.nComps);

% Pre allocation for GLMs
out.StudyGLMs       		= cell(out.nChans,1);
out.StudyGLMSChanCompTVal 	= zeros(out.nChans,out.nComps);
out.StudyGLMsChanRsquared   = zeros(out.nChans,1);

out.TestGLMs        		= cell(out.nChans,1);
out.TestGLMSChanCompTVal 	= zeros(out.nChans,out.nComps);
out.TestGLMsChanRsquared    = zeros(out.nChans,1);

for ss=1:out.nSubjs
    subjChans = find(data.subjChans==ss);
	% re-order to channels, trials , bands , time bins
    x = permute(data.BinERPs{ss}(:,:,:,data.AnalysisBins),[2 3 1 4]);
    % concatenate bands and time bins
    x = x(:,:,:);
    [d1,d2,d3] = size(x); % d1 -> n subj channels, d2 -> # of trials, d3 -> bands*timebins
    rts1 = data.studyRTs{ss}(data.trials{ss});
    rts2 = data.testRTs{ss}(data.trials{ss});

    out.Comps{ss} = zeros(d1,d2,out.nComps);

        
    % for subject channels
    for c=1:d1
    	% get the trial by band*bin matrix for each channel.
        xx = squeeze(x(c,:,:));

        % PCA
        [C,S,E]= ppca(xx',out.nComps);  
        out.Comps{ss}(c,:,:) = C;
        out.Projections(subjChans(c),:,:) = S;
        out.VarExp(subjChans(c),:) = E;
        
        % Correlations
        out.CorrStudyRTs(subjChans(c),:) = corr(C,atan(rts1));
        out.CorrTestRTs(subjChans(c),:)  = corr(C,atan(rts2));
        
        % GLMSs
		out.StudyGLMs{subjChans(c)} = fitglm(C,atan(rts1));
		out.StudyGLMSChanCompTVal(subjChans(c),:) ...
			= out.StudyGLMs{subjChans(c)}.Coefficients.tStat(2:out.nComps+1);
		out.StudyGLMsChanRsquared(subjChans(c))...
			= out.StudyGLMs{subjChans(c)}.Rsquared.Ordinary;

		out.TestGLMs{subjChans(c)} = fitglm(C,atan(rts2));
		out.TestGLMSChanCompTVal(subjChans(c),:) ...
			= out.TestGLMs{subjChans(c)}.Coefficients.tStat(2:out.nComps+1);
		out.TestGLMsChanRsquared(subjChans(c)) ...
			= out.TestGLMs{subjChans(c)}.Rsquared.Ordinary;
    end
end