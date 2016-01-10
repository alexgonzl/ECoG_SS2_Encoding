function [SpCorr,SpPval] = corrDataToRTs(data,RTs)

if size(RTs,2)>2
    error('RTs not in right format');
end
if size(RTs,1)~=size(data,1);
    error('data sizes do not match')
end
notValidRTs = isnan(RTs);
RTs(notValidRTs) = [];
data(notValidRTs,:)=[];

% convert RTs to log scale
RTs = log10(RTs);


%% get main values
[SpCorr,SpPval] = corr(data,RTs,'type','spearman');


%% permutation dist

% B = 1000;
% N = numel(RTs);   % number of samples
% P = size(data,2); % number of instances
% cB = zeros(B,P);
% for bb = 1:B    
%     randOrd = randperm(N);
%     cB(bb,:) = corr(data(randOrd,:),RTs,'type','spearman');
% end
% cB = sort(cB);
% 
% out.PermPval = zeros(P,1);
% for pp = 1:P
%     y=cB(:,pp);
%     x=out.SpCorr(pp);    
%     
%     % two sided test
%     if x>0
%         ngreaterthan = sum(x>y);
%         out.PermPval(pp) = (1-ngreaterthan/B)*2;
%     else
%         nlessthan    = sum(x<y);
%         out.PermPval(pp) = (1-nlessthan/B)*2;
%     end
%     % default to lowest possible # based # of permutations
%     if out.PermPval(pp)==0
%         out.PermPval(pp) = 2/B;
%     end
% end
