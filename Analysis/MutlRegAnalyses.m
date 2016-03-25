function out = MultRegAnalyses(data,opts)
% function performs multiple regression on a subject level 
% trying to predict individuals subjects reactions times
%
% 	opts:
%		regType   	-> {'lasso', 'ridge'}
% 		multiChan 	-> boolean. if true, uses all available channels for each subejct.
%		KFold 	  	-> number of XVal folds.
% 		RTset 		-> {'studyRTs','testRTs'}




% % use multiple regularized regression for testCorrs
% data.RidgeTestRTs.nXval = 10;
% data.RidgeTestRTs.Y     = cell(nSubjs,1);
% data.RidgeTestRTs.Yhat  = cell(nSubjs,1);
% data.RidgeTestRTs.resid  = cell(nSubjs,1);
% data.RidgeTestRTs.corr   = zeros(nSubjs,1);
% data.RidgeTestRTs.TestIDx = cell(nSubjs,1);
% data.RidgeTestRTs.GenModel = zeros(data.nChans,data.nBands, data.nBins);

% data.LassoTestRTs.nXval = 10;
% data.LassoTestRTs.B = cell(nSubjs,1);
% data.LassoTestRTs.STATS = cell(nSubjs,1);
% data.LassoTestRTs.GenModel = zeros(data.nChans,data.nBands,data.nBins);
% data.LassoTestRTs.Yhat  = cell(nSubjs,1);
% data.LassoTestRTs.Y     = cell(nSubjs,1);
% data.LassoTestRTs.corr  = zeros(nSubjs,1);
% for ss = 1:nSubjs
%     N = sum(data.trials{ss}); 
%     data.nTrialsInAnalysis(ss)=N;
%     data.LassoTestRTs.Yhat{ss} = zeros(N,1);
%     data.LassoTestRTs.TestIDx{ss} = crossvalind('Kfold', N, data.RidgeTestRTs.nXval);    

% %     data.RidgeTestRTs.Yhat{ss} = zeros(N,1);
% %     data.RidgeTestRTs.TestIDx{ss} = crossvalind('Kfold', N, data.RidgeTestRTs.nXval);    
% %     
%     X=data.BinERPs{ss}(:,:);
%     Y=data.testRTs{ss}(data.trials{ss}); 
% %     data.RidgeTestRTs.Y{ss} = Y;
%     data.LassoTestRTs.Y{ss} = Y;
%     for kk = 1:data.LassoTestRTs.nXval 
%          IDx=data.LassoTestRTs.TestIDx{ss}~=kk;
%          [B,STATS] = lasso(X(IDx,:),Y(IDx),'alpha',numel(Y)/numel(X),'CV',10,'DFMax',numel(Y));
%          Bo = B(:,STATS.IndexMinMSE);
%          data.LassoTestRTs.Yhat{ss}(~IDx) = X(~IDx,:)*Bo;
%     end
%     data.LassoTestRTs.corr(ss) = corr(data.LassoTestRTs.Yhat{ss},data.LassoTestRTs.Y{ss});
    
% %     temp = ridge(Y,X,numel(X)/numel(Y),0);
% %     temp(1) = [];
% %     data.RidgeTestRTs.GenModel(data.subjChans==ss,:,:) = reshape(temp,[data.nSubjLPCchans(ss),data.nBands,data.nBins]);
% %     for kk = 1:data.RidgeTestRTs.nXval        
% %         IDx=data.RidgeTestRTs.TestIDx{ss}~=kk;
% %         b0 = ridge(Y(IDx),X(IDx,:),numel(X)/numel(Y),0);
% %         data.RidgeTestRTs.Yhat{ss}(~IDx) = [ones(sum(~IDx),1) ...
% %             X(~IDx,:)]*b0;       
% %     end
% %    data.RidgeTestRTs.resid{ss}= Y-data.RidgeTestRTs.Yhat{ss};
% %    data.RidgeTestRTs.corr(ss) = corr(Y,data.RidgeTestRTs.Yhat{ss});
%    fprintf(' xval completed for subj %i',ss);
   
% %     % Lasso 
% %     [B,STATS] = lasso(X,Y,'alpha',numel(Y)/numel(X),'CV',10,'DFMax',numel(Y));      
% % 
% %     data.LassoTestRTs.B{ss} = B;
% %     data.LassoTestRTs.STATS{ss} = STATS;
% %     data.LassoTestRTs.GenModel(data.subjChans==ss,:,:) = ...
% %         reshape(B(:,STATS.IndexMinMSE),[data.nSubjLPCchans(ss),data.nBands,data.nBins]);    
% %     data.LassoTestRTs.Yhat{ss} = X*B(:,STATS.IndexMinMSE);
% %     data.LassoTestRTs.corr(ss) = corr(data.LassoTestRTs.Yhat{ss},data.LassoTestRTs.Y{ss});
% end

% %%
% % use multiple regularized regression for study RTs
% data.RidgeStudyRTs.nXval = 10;
% data.RidgeStudyRTs.Y     = cell(nSubjs,1);
% data.RidgeStudyRTs.Yhat  = cell(nSubjs,1);
% data.RidgeStudyRTs.resid  = cell(nSubjs,1);
% data.RidgeStudyRTs.corr   = zeros(nSubjs,1);
% data.RidgeStudyRTs.TestIDx = cell(nSubjs,1);
% data.RidgeStudyRTs.GenModel = zeros(data.nChans,data.nBands,data.nBins);

% data.LassoStudyRTs.nXval = 10;
% data.LassoStudyRTs.B = cell(nSubjs,1);
% data.LassoStudyRTs.STATS = cell(nSubjs,1);
% data.LassoStudyRTs.GenModel = zeros(data.nChans,data.nBands,data.nBins);
% data.LassoStudyRTs.Yhat  = cell(nSubjs,1);
% data.LassoStudyRTs.Y  = cell(nSubjs,1);
% data.LassoStudyRTs.corr  = zeros(nSubjs,1);

% for ss = 1:nSubjs
%     %N = sum(data.conds{ss,10}); 
% %     data.RidgeStudyRTs.Yhat{ss} = zeros(N,1);
% %     data.RidgeStudyRTs.TestIDx{ss} = crossvalind('Kfold', N, data.RidgeStudyRTs.nXval);    
% %     
%     X=data.BinERPs{ss}(:,:);
%     Y=data.studyRTs{ss}(data.trials{ss});
%     data.LassoStudyRTs.Y{ss} = Y;
    
%     N = sum(data.trials{ss}); 
%     data.nTrialsInAnalysis(ss)=N;
%     data.LassoStudyRTs.Yhat{ss} = zeros(N,1);
%     data.LassoStudyRTs.TestIDx{ss} = crossvalind('Kfold', N, data.LassoStudyRTs.nXval);    

% %     data.RidgeTestRTs.Yhat{ss} = zeros(N,1);
% %     data.RidgeTestRTs.TestIDx{ss} = crossvalind('Kfold', N, data.RidgeTestRTs.nXval);    
% %        
% %     data.RidgeTestRTs.Y{ss} = Y;
%     data.LassoTestRTs.Y{ss} = Y;
%     for kk = 1:data.LassoStudyRTs.nXval 
%          IDx=data.LassoStudyRTs.TestIDx{ss}~=kk;
%          [B,STATS] = lasso(X(IDx,:),Y(IDx),'alpha',numel(Y)/numel(X),'CV',10,'DFMax',numel(Y));
%          Bo = B(:,STATS.IndexMinMSE);
%          data.LassoStudyRTs.Yhat{ss}(~IDx) = X(~IDx,:)*Bo;
%     end
%     data.LassoStudyRTs.corr(ss) = corr(data.LassoStudyRTs.Yhat{ss},data.LassoStudyRTs.Y{ss});
    
% %     temp = ridge(Y,X,numel(X)/numel(Y),0);
% %     data.RidgeStudyRTs.Y{ss} = Y;
% %     temp = ridge(Y,X,numel(X)/numel(Y),0);
% %     temp(1) = [];
% %     data.RidgeStudyRTs.GenModel(data.subjChans==ss,:,:) = ...
% %         reshape(temp,[data.nSubjLPCchans(ss),data.nBands,data.nBins]);
% %     for kk = 1:data.RidgeStudyRTs.nXval
% %         IDx=data.RidgeStudyRTs.TestIDx{ss}~=kk;
% %         b0 = ridge(Y(IDx),X(IDx,:),numel(X)/numel(Y),0);
% %         data.RidgeStudyRTs.Yhat{ss}(~IDx) = [ones(sum(~IDx),1) ...
% %             X(~IDx,:)]*b0;
% %     end
% %     data.RidgeStudyRTs.resid{ss}= Y-data.RidgeStudyRTs.Yhat{ss};
% %     data.RidgeStudyRTs.corr(ss) = corr(Y,data.RidgeStudyRTs.Yhat{ss});
%     fprintf(' xval completed for subj %i',ss);
    
% %     % Lasso 
% %     [B,STATS] = lasso(X,Y,'alpha',numel(Y)/numel(X),'CV',10,'DFMax',numel(Y));      
% %     data.LassoStudyRTs.B{ss} = B;
% %     data.LassoStudyRTs.STATS{ss} = STATS;
% %     data.LassoStudyRTs.GenModel(data.subjChans==ss,:,:) = ...
% %         reshape(B(:,STATS.IndexMinMSE),[data.nSubjLPCchans(ss),data.nBands,data.nBins]);    
% %     data.LassoStudyRTs.Yhat{ss} = X*B(:,STATS.IndexMinMSE);
% %     data.LassoStudyRTs.corr(ss) = corr(data.LassoStudyRTs.Yhat{ss},data.LassoStudyRTs.Y{ss});
% end



