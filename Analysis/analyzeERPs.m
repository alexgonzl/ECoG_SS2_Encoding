%function data=analyzeERPs(data);

% first get a split of RT quantiles
data.RTquantiles = [0.4 0.6];
remTrials = data.behavior.subRem==1; % correctly remember trials
data.testRTsQuants = quantile(data.behavior.testRTs(remTrials),data.RTquantiles);

data.subFastRem     = remTrials & (data.behavior.testRTs <= data.testRTsQuants(1));
data.subSlowRem     = remTrials & (data.behavior.testRTs >= data.testRTsQuants(2));


%%
IPSch = data.ROIid==1;
SPLch = data.ROIid==2;

X   = data.BinZStat(IPSch,:);
Y   = data.BinZStat(SPLch,:);

[~,pIPS,~,tIPS]=ttest(X);
[~,pSPL,~,tSPL]=ttest(Y);
[~,pInt,~,tInt]=ttest2(X,Y);