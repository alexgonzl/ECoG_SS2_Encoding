function out = getTwoSampStats(x,y)
% produces several two sample statistics
% 
% the different statistics are represented on rows in the following order:
% 1) t statstics (same mean)
% 2) Kolmogorov-Smirnov (same distribution)
% 3) Wilcoxon (same medians)
%
% the firs column indicates the p value, and the second indicates the
% resulting statistic


x(isnan(x)) = [];
y(isnan(y)) = [];

out = zeros(3,2);
[~,temp1,~,temp2]=ttest2(x,y);
out(1,1) = temp1;
out(1,2) = temp2.tstat;
[~,out(2,1),out(2,2)] = kstest2(x,y);
[out(3,1), out(3,2)] = ranksum2(x,y);