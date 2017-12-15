function PrePostActRT_RegionPairs_Plots(opts)

fileName = ['allMBAnalysis_2preStimsublogPowernonLPCch'];
load([opts.dataPath opts.lock '/' fileName '.mat'])
% fileName = ['allLME_ModelAnalysis' opts.lock 'sublogPowernonLPCch'];
% load([opts.dataPath opts.lock '/' fileName '.mat'])
% data=out;

inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';
%fP = '~/Google Drive/Research/ECoG_SS2e/Plots/Behavior/';
fP = opts.savePath;

ROIcolors = [240 35 17; 2 93 140;122 179 23]/255;

% summary of data.
post_str = '_1500ms_poststim';

% Display of Results:
fprintf('\n')
disp('Study RT to Pre/Post Activity to models')
disp(data.Chan2ChanPrePostAct2RT.studyModel.Formula)
disp(data.Chan2ChanPrePostAct2RT.studyANOVA)

fprintf('\n')
disp('Test RT to Pre/Post Activity to models')
disp(data.Chan2ChanPrePostAct2RT.testModel.Formula)
disp(data.Chan2ChanPrePostAct2RT.testANOVA)


% var names
conds_s = {'pre','post','pre->post'}; nConds = numel(conds_s);
roi_pairs  = {'IPS-IPS','IPS-SPL','IPS-AG','SPL-IPS','SPL-SPL',...
    'SPL-AG','AG-IPS','AG-SPL','AG-AG'};
roi_pairs2  = {'IPS','SPL','AG','IPS','SPL','AG','IPS','SPL','AG'};
rois ={'pre-IPS', 'pre-SPL', 'pre-AG'};
nROIs = numel(rois);


%% activity correlation plot

close all;
mainEff = {'conds', 'rois', 'bands'};
%xTickLabels = {conds_s,rois_s,bands_s2};
%varNames = {conds_s,rois_s,bands_s};

Colors{1} = brewermap(4,'Reds');
Colors{1}(1,:) = [];
Colors{2} = brewermap(4,'Blues');
Colors{2}(1,:) = [];
Colors{3} = brewermap(4,'Greens');
Colors{3}(1,:) = [];

RTs_s = {'study','test'};
conds = {'PrePostC'};
nConds = numel(conds);
xTickLabelSet = [];
xTickLabelSet{1} =  roi_pairs2;
xTickLabelSet{2} =  rois;

Tb = data.Chan2ChanPrePostAct2RT.table;

if 1
    for ii = 1:nConds % study/test
        X = cell(nROIs);
        cnt = 1;
        for r1 =1:nROIs % interactions
            for r2 = 1:nROIs
                X{r1,r2} = tableSlice(Tb,conds{ii},'ROI_Pairs', ...
                    roi_pairs{cnt});
                    cnt=cnt+1;
            end
        end
        fN =  strcat('Chan2ChanRegionEffects_',conds{ii},'_',post_str, '_2');
        groupBarPlot(X,Colors,xTickLabelSet,strcat(fP,fN))
    end
end


%% interactions with RTs
conds = {'T3','S3'};
prepostSign = {'pos','neg'};
Colors{1} = brewermap(4,'Reds');
Colors{1}(1,:) = [];
Colors{2} = brewermap(4,'Blues');
Colors{2}(1,:) = [];
Colors{3} = brewermap(4,'Greens');
Colors{3}(1,:) = [];

% table excluding low corr channel pairs.
Tb = data.Chan2ChanPrePostAct2RT.table2;


for ii = 1:2
for jj = 1:2
    X=cell(nROIs);
    cnt = 1;
for r1=1:nROIs
    for r2 =1:nROIs 
        X{r1,r2}=tableSlice2(Tb, conds{ii},'ROI_Pairs', roi_pairs{cnt} ,'PrePostC2',prepostSign{jj});
        cnt=cnt+1;
    end
end
    fN =  strcat('Chan2ChanRegionEffects_',prepostSign{jj},'_',conds{ii},'_',post_str);
        groupBarPlot(X,Colors,xTickLabelSet,strcat(fP,fN))
end
end

 end                       
%%
function t = ttest_T(dat)
[~,~,~,t] = ttest(dat);
t = t.tstat;
end

function p = ttest_P(dat)
[~,p] = ttest(dat);
end

function datSlice = tableSlice(dataTable,outVarName,colName,varName)
idx = strcmp(dataTable.(colName),varName);
datSlice = dataTable.(outVarName)(idx);
end

function datSlice = tableChanSlice(dataTable,outVarName,colName,varName)
idx = strcmp(dataTable.(colName),varName);
c = string(dataTable.chans(idx));
x = dataTable.(outVarName)(idx);
datSlice = grpstats(x,c,'mean');

end

function datSlice = tableSlice2(dataTable,outVarName,colName1,varName1,colName2,varName2)
idx = strcmp(dataTable.(colName1),varName1);
idx = idx & strcmp(dataTable.(colName2),varName2);
datSlice = dataTable.(outVarName)(idx);
end

function datSlice = tableChanSlice2(dataTable,outVarName,colName1,varName1,colName2,varName2)
idx = strcmp(dataTable.(colName1),varName1);
idx = idx & strcmp(dataTable.(colName2),varName2);
c = string(dataTable.chans(idx));
x = dataTable.(outVarName)(idx);
datSlice = grpstats(x,c,'mean');

end

function datSlice = tableChanSlice3(dataTable,outVarName,colName1,varName1,colName2,varName2,colName3,varName3)
idx = strcmp(dataTable.(colName1),varName1);
idx = idx & strcmp(dataTable.(colName2),varName2);
idx = idx & strcmp(dataTable.(colName3),varName3);
c = string(dataTable.chans(idx));
x = dataTable.(outVarName)(idx);
datSlice = grpstats(x,c,'mean');

end


function simpBarPlot(M,colors,xticklabel,fN)

N = numel(M);
opts =[];
opts.aspectRatio = [100*N 400];
opts.xTick = 1:N;
opts.xTickLabel = xticklabel;
opts.yTicks = 'def';
opts.grid = 1;
opts.yLabel = 'Across Channel Means (Ts)';
opts.colors = colors;

pp=cellTStatWrapper(M);
[~,a] = fdr_bh(pp.uniVariateP);
ids=pp.uniVariateP<=a;
opts.sigUniMarks = find(ids);
opts.sigUniLevel = pp.uniVariateP(ids);

barPlotWithErrors(M,opts);
print(gcf,'-dpdf', strcat(fN,'b'))
scatterColPlot(M,opts);
print(gcf,'-dpdf', strcat(fN,'s'))
end


function groupBarPlot(X,colors,xticklabel,fN)

[N,M] = size(X);

opts =[];
opts.aspectRatio = [90*N*M 400];
opts.xTick = 'def';
opts.xTickLabel = xticklabel;
opts.yTicks = 'def';
opts.grid = 1;
opts.yLabel = 'Across Channel Pairs Means (Rs)';
opts.colors = colors;
opts.scatter = 's';

p=cellfun(@ttest_P,X);
[~,a] = fdr_bh(p(:));

opts.sigUniMarks = p<=a;
opts.sigUniLevel = p;

ha=figure(); clf;hold on;
set(gcf,'paperpositionmode','auto','position',[100 100 opts.aspectRatio], 'papersize',[15 6])
barPlotWithErrors_grp(X,opts);
print(gcf,'-dpdf',strcat(fN,'b',opts.scatter))


% ha=figure(); clf;hold on;
% set(gcf,'paperpositionmode','auto','position',[100 100 opts.aspectRatio], 'papersize',[20 10])
% scatterColPlot_grp(X,opts);
% print(gcf,'-dpdf', strcat(fN,'s'))
end

function saveFig(figNum, fP, fN)
% save plot
inkscape ='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';

try
    cPath = pwd;
    cd(fP)
    addpath(cPath)
    addpath([cPath '/Plotting/'])
    
    print(figNum,'-dsvg',fN)
    eval(['!' inkscape ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
    eval(['!rm ' fN '.svg'])
    
    cd(cPath)
catch
    cd(cPath)
end
end
