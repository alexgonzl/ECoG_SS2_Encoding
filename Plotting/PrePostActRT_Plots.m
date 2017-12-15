function PrePostActRT_Plots(opts)
%fileName = ['allMBAnalysis_2' opts.lock 'sublogPowernonLPCch'];
%load([opts.dataPath opts.lock '/' fileName '.mat'])
fileName = ['allLME_ModelAnalysis' opts.lock 'sublogPowernonLPCch'];
load([opts.dataPath opts.lock '/' fileName '.mat'])
data=out;

inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';
%fP = '~/Google Drive/Research/ECoG_SS2e/Plots/Behavior/';
fP = opts.savePath;

ROIcolors = [240 35 17; 2 93 140;122 179 23]/255;

% summary of data.
if opts.poststim500==1
    post_str = '_500ms_poststim';
    Tb = data.PrePostActRT.datTable2;
    
    % Display of Results:
    disp('P val threshold for MC ANOVAs for alpha=0.05')
    disp(data.PrePostActRT.AcrossANOVAs_q2)
    fprintf('\n')
    disp('Study RT to Pre/Post Activity to models')
    disp(data.PrePostActRT.studyModel2.Formula)
    disp(data.PrePostActRT.studyANOVA_FDR2)
    disp(data.PrePostActRT.studyModel2.ModelCriterion)
    
    fprintf('\n')
    disp('Test RT to Pre/Post Activity to models')
    disp(data.PrePostActRT.testModel2.Formula)
    disp(data.PrePostActRT.testANOVA_FDR2)
    disp(data.PrePostActRT.testModel2.ModelCriterion)
    
else
     post_str = '_1500ms_poststim';
    Tb = data.PrePostActRT.datTableX;
    
    % Display of Results:
    disp('P val threshold for MC ANOVAs for alpha=0.05')
    disp(data.PrePostActRT.AcrossANOVAs_qX)
    fprintf('\n')
    disp('Study RT to Pre/Post Activity to models')
    disp(data.PrePostActRT.studyModelX.Formula)
    disp(data.PrePostActRT.studyANOVA_FDRX)
    disp(data.PrePostActRT.studyModelX.ModelCriterion)
    
    fprintf('\n')
    disp('Test RT to Pre/Post Activity to models')
    disp(data.PrePostActRT.testModelX.Formula)
    disp(data.PrePostActRT.testANOVA_FDRX)
    disp(data.PrePostActRT.testModelX.ModelCriterion)
   
end

% var names
conds_s = {'pre','post'}; nConds = numel(conds_s);
rois_s  = {'IPS','SPL','AG'};  nROIs = numel(rois_s);
bands_s = {'delta','theta','alpha','beta','lgam','hgam'}; nBands = numel(bands_s);
bands_s2 = {'\delta','\theta','\alpha','\beta','l-\gamma','h-\gamma'};



%% Main Effects Plots
close all;
mainEff = {'bands', 'conds', 'rois'};
xTickLabels = {bands_s2,conds_s,rois_s};
varNames = {bands_s,conds_s,rois_s};

RTs_s = {'study','test'};
% color Sets:
close all
Colors{1} = brewermap(6,'spectral');
Colors{2} = brewermap(3,'Pastel2');
Colors{3} = ROIcolors;

if 0
    for ii = 1:2 % study/test
        for jj =1:3 % main Eff
            N = numel(varNames{jj});
            X = cell(N,1); Y = cell(N,1);
            for kk = 1:N
                X{kk} = tableSlice(Tb,RTs_s{ii},mainEff{jj},varNames{jj}{kk});
                Y{kk} = tableChanSlice(Tb,RTs_s{ii},mainEff{jj},varNames{jj}{kk});
            end
            
            fN = strcat('MainEffects_',RTs_s{ii},'_',mainEff{jj}, post_str, '_1');
            simpBarPlot(X,Colors{jj},xTickLabels{jj},strcat(fP,fN));
            fN = strcat('MainEffects_',RTs_s{ii},'_',mainEff{jj},post_str, '_2');
            simpBarPlot(Y,Colors{jj},xTickLabels{jj},strcat(fP,fN));
            
        end
    end
end
%% Two level interaction plots

close all;
mainEff = {'conds', 'rois', 'bands'};
xTickLabels = {conds_s,rois_s,bands_s2};
varNames = {conds_s,rois_s,bands_s};

Colors{1} = brewermap(3,'Pastel2');
Colors{2} = ROIcolors;
Colors{3} = brewermap(6,'spectral');

RTs_s = {'study','test'};

if 0
    for ii = 1:2 % study/test
        for jj =1:2 % interactions
            N = numel(varNames{jj});
            for kk = (jj+1):3
                M = numel(varNames{kk});
                X = cell(N,M);
                Y = cell(N,M);
                ColorSet = cell(N,1);
                for j = 1:N
                    for k = 1:M
                        X{j,k} = tableSlice2(Tb,RTs_s{ii},mainEff{jj}, ...
                            varNames{jj}{j},mainEff{kk},varNames{kk}{k});
                        Y{j,k} = tableChanSlice2(Tb,RTs_s{ii},mainEff{jj}, ...
                            varNames{jj}{j},mainEff{kk},varNames{kk}{k});
                        
                    end
                    ColorSet{j} = Colors{kk};
                end
                xTickLabelSet{1} =  xTickLabels{kk};
                xTickLabelSet{2} =  xTickLabels{jj};
                
                fN =  strcat('InterEffects_',RTs_s{ii},'_', mainEff{jj}, '_',  mainEff{kk},post_str, '_1');
                groupBarPlot(X,ColorSet,xTickLabelSet,strcat(fP,fN))
                fN =  strcat('InterEffects_',RTs_s{ii},'_', mainEff{jj}, '_',  mainEff{kk},post_str, '_2');
                groupBarPlot(Y,ColorSet,xTickLabelSet,strcat(fP,fN))
            end
        end
    end
end

%% Three way interaction plot

close all;
mainEff = {'conds', 'rois', 'bands'};
xTickLabels = {conds_s,rois_s,bands_s2};
varNames = {conds_s,rois_s,bands_s};

Colors{1} = brewermap(3,'Pastel2');
Colors{2} = ROIcolors;
%Colors{3} = brewermap(6,'spectral');
Colors{3} = brewermap(12,'BuPu');
Colors{3}(1:6,:)=[];

RTs_s = {'study','test'};
N = numel(conds_s);
M = numel(rois_s);
P = numel(bands_s);

nConds=2;
ColorSet = cell(M,1);
for j = 1:M
    ColorSet{j} = Colors{3};
end

nlm = cell(2,2,3);
nlm2 = cell(2,2,3);
nlm3 = cell(2,2,3);
nlm4 = cell(2,2,3);

nlmT = zeros(3,2,2,3);
nChans = [59,49,26];
modelfun = @(b,x)(b(1) + b(2)*x +b(3)*(x.^2) );
modelfun2 = @(b,x)(b(1)  + b(2)*(x.^2) );
modelfun3 = @(b,x)(b(1)  + b(2)*x  );
modelfun4 = @(b,x)(b(1) + x );
nlmBIC = zeros(4,2,2,3);

for ss = 1:2
    X = cell(N,M,P);
    Y = cell(N,M,P);
    for i= 1:N
        for j = 1:M
            for k = 1:P
                Y{i,j,k} = tableChanSlice3(Tb,RTs_s{ss},'conds', ...
                    conds_s{i},'rois',rois_s{j},'bands',bands_s{k});                
            end
            y=[Y{i,j,1:P}]; y=y(:);
            x=repmat(1:P,[nChans(j),1]);x=x(:);

            nlm{ss,i,j}=fitnlm(x,y,modelfun,[0 1 1]);                        
            nlm2{ss,i,j}    =fitnlm(x,y,modelfun2,[0 1]);
            nlm3{ss,i,j}    =fitnlm(x,y,modelfun3,[0 1]);
            nlm4{ss,i,j}    =fitnlm(x,y,modelfun4,[0]);
            nlmT(:,ss,i,j)  = nlm{ss,i,j}.Coefficients.tStat;
            nlmBIC(1,ss,i,j)  = nlm{ss,i,j}.ModelCriterion.BIC;
            nlmBIC(2,ss,i,j)  = nlm2{ss,i,j}.ModelCriterion.BIC;
            nlmBIC(3,ss,i,j)  = nlm3{ss,i,j}.ModelCriterion.BIC;
            nlmBIC(4,ss,i,j)  = nlm4{ss,i,j}.ModelCriterion.BIC;
            
        end
    end
    
    xTickLabelSet{1} =  bands_s2;
    xTickLabelSet{2} =  rois_s;
    for ii = 1:nConds
        fN =  strcat('TripInterEffects_',RTs_s{ss},'_', conds_s{ii},post_str, '_rois_bands');
        X = squeeze(Y(ii,:,:));
%        groupBarPlot(X,ColorSet,xTickLabelSet,strcat(fP,fN))
    end
    
    
end

x=Y;
end

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
opts.aspectRatio = [70*N*M 400];
opts.xTick = 'def';
opts.xTickLabel = xticklabel;
opts.yTicks = 'def';
opts.grid = 1;
opts.yLabel = 'Across Channel Means (Ts)';
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
