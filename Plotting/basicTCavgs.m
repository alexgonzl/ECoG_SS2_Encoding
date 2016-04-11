
load('/Volumes/ECoG_SS2/SS2/SS2e/Results/group/Spectral_Data/allERSPshgamGroupstimLocksublogPowernonLPCCh.mat')
rRTs = log10(cell2mat(data.testRTs'));
IPSch = data.hemChan==1 & data.ROIid==1;
SPLch = data.hemChan==1 & data.ROIid==2;
%% Activity in IPS split by fast/slow

X = [];
X{1}=data.testmCond1(IPSch,:);
X{2}=data.testmCond2(IPSch,:);
plotNTraces(X,data.trialTime,'oc','loess',0.15)

%% Activity in SPL split by fast/slow
X = [];
X{1}=data.testmCond1(SPLch,:);
X{2}=data.testmCond2(SPLch,:);
plotNTraces(X,data.trialTime,'oc','loess',0.15)

%% Activity in IPS split by abs conc

X = [];
X{1}=data.studymCond1(IPSch,:);
X{2}=data.studymCond2(IPSch,:);
plotNTraces(X,data.trialTime,'oc','loess',0.15)

%% Activity in SPL split by abs conc

X = [];
X{1}=data.studymCond1(SPLch,:);
X{2}=data.studymCond2(SPLch,:);
plotNTraces(X,data.trialTime,'oc','loess',0.15)


%% IPS and SPL RT correlation to activity

X = [];
X{1}=data.BinRT_RemTrialsZCorr(IPSch,:);
X{2}=data.BinRT_RemTrialsZCorr(SPLch,:);

plotNTraces(X,t,'rb')


%% IPS: activity to test RT correlation split by abs concrete
X = [];
X{1}=data.BinStudyZcCond1(IPSch,:);
X{2}=data.BinStudyZcCond2(IPSch,:);

plotNTraces(X,t,'oc')

%% SPL: activity to test RT correlation split by abs concrete
X = [];
X{1}=data.BinStudyZcCond1(SPLch,:);
X{2}=data.BinStudyZcCond2(SPLch,:);

plotNTraces(X,t,'oc')

