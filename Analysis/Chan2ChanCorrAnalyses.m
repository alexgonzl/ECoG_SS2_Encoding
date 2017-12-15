
function data = Chan2ChanCorrAnalyses(data)
%%
nSubjs = 8;
preBins = 3;
postBins = 4:6;

nP = 2; % predictions are to study/test RTs
nConds = 3;
ROIs = {'IPS','SPL','AG'};
nROIs=3;

T = nan(nROIs,nROIs,nConds,nP,nSubjs);
P = nan(nROIs,nROIs,nConds,nP,nSubjs);

Study = cell(nROIs,1);
Test = cell(nROIs,1);

R_Vec = [];
Hem_Vec = [];

ROI_Combs = [];
Subj_Vec  = [];
for ss = 1:nSubjs
    
    subjROIchans = data.ROIid(data.subjChans==ss);
    trials = data.trials{ss};
    
    RTs = [];
    RTs{1} = log10(data.studyRTs{ss}(trials));
    RTs{2} = log10(data.testRTs{ss}(trials));
    
    nTrials = sum(trials);
    nChans = data.nSubjLPCchans(ss);
    % conds: pre post pre to post (#3)
    % predicting: studyRTs/testRTs (#2)
    bivariateChan_Act_To_RTCorr = nan(nChans,nChans,nConds,nP,nSubjs);
    
    x=squeeze(mean(data.BinERPs2{ss}(:,:,:,preBins),4));
    y=squeeze(mean(data.BinERPs2{ss}(:,:,:,postBins),4));
    
    n2 = nChans^2;
    n2_2 = (n2-nChans)/2;
    mask = tril(true(nChans),-1);
    
    C = cell(3,1);
    C{1} = zeros(nTrials,n2_2);
    C{2} = zeros(nTrials,n2_2);
    C{3} = zeros(nTrials,n2);
    
    M = cell(nChans,nChans);
    for ii=1:nChans
        for jj=1:nChans
            M{ii,jj} = strcat(ROIs{subjROIchans(ii)},'-',ROIs{subjROIchans(jj)});
        end
    end
    for ii = 1:nTrials
        xx=x(:,:,ii);
        yy=y(:,:,ii);
        c1 = atanh(corr(xx));
        c2 = atanh(corr(yy));
        c3 = atanh(corr(xx,yy));
        C{1}(ii,:) = c1(mask)';
        C{2}(ii,:) = c2(mask)';
        C{3}(ii,:) = c3(:)';
    end
    
    for pp =1:nP % study and test RTs
        for cc = 1:nConds % number of conditions (pre/post pre to post)
            if cc==3
                 temp=reshape(corr(C{cc},RTs{pp},'type','spearman'),[nChans,nChans]);
                 if pp==1
                     R_Vec       = [R_Vec;repmat(subjROIchans,[nChans,1])];
                    ROI_Combs   = [ROI_Combs;M(:)];
                    Subj_Vec    = [Subj_Vec;ss*ones(nChans^2,1)];
                    Hem_Vec     = [Hem_Vec; repmat(data.hemChan(data.subjChans==ss),[nChans,1])];
                 end
            else
                temp = nan(nChans);
                temp(mask) =corr(C{cc},RTs{pp},'type','spearman');
            end
            if pp==1
                Study{cc}   = [Study{cc};temp(:)];                
             else
                 Test{cc} = [Test{cc}; temp(:)];
             end
            bivariateChan_Act_To_RTCorr(:,:,cc,pp,ss) = temp;
            
        end
    end

    for pp = 1:nP
        for cc = 1:nConds
            for r1=1:nROIs
                if any(subjROIchans==r1)
                    for r2=1:nROIs
                        if any(subjROIchans==r2)
                            x = bivariateChan_Act_To_RTCorr(subjROIchans==r1,subjROIchans==r2,cc,pp,ss);
                            if numel(x)>1
                                [~,p,~,t] = ttest((x(:)));
                                
                                T(r1,r2,cc,pp,ss) = t.tstat;
                                P(r1,r2,cc,pp,ss) = p;
                            end
                        end
                    end
                end
            end
        end
    end
    
end
tbl = table(Study{1},Study{2},Study{3},Test{1},Test{2},Test{3},...
    ROI_Combs, categorical(Hem_Vec), categorical(Subj_Vec),'VariableNames',...
    {'S1','S2','S3','T1','T2','T3','ROI_Pairs','Hem','Subjs'});
%%
% TS = nan(nROIs,nROIs,nConds,nP);
% PS = nan(nROIs,nROIs,nConds,nP);
% LeftChans = unique(data.subjChans(data.hemChan==1));
% for pp = 1:nP
%     for cc = 1:nConds
%         for r1=1:nROIs
%             for r2=1:nROIs
%                 x = squeeze(T(r1,r2,cc,pp,:));
%                 %x(LeftChans) = -x(LeftChans);
%                 
%                 [~,p,~,t] = ttest(x);
% 
%                 TS(r1,r2,cc,pp) = t.tstat;
%                 PS(r1,r2,cc,pp) = p;
%             end
%         end
%     end
% end
end
%%
function t = ttest_T(dat)
[~,~,~,t] = ttest(dat);
t = t.tstat;
end

function p = ttest_P(dat)
[~,p] = ttest(dat);
end


