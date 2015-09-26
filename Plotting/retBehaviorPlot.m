% script for plotting behavioral performance at retrieval

X = [0.91	0.87	0.92	0.96	0.88	0.91
0.14	0.12	0.01	0.04	0.24	0.11
2.41	2.28	3.66	3.54	1.89	2.76
0.97	0.99	0.97	0.88	0.53	0.87
0.77	0.96	1.17	1.1	1.14	1.03
1.28	1.01	1.07	1.6	1.46	1.28];


%%
shapes = 'odsph';
figure(1); clf; set(gcf, 'position', [100 100 1000 300],'paperpositionmode','auto'); 
sD = 200;

%
% Hit Rate
axes('position',[0.05 0.1 0.15 0.85]); hold on;
xlim([0 1]); ylim([-0.1 1.1])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,X(1,ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD)
end
plot([0.2 0.8],[1 1]*X(1,nSubjs+1),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 0.5 1],'linewidth',2); 
xlabel(' Hit Rate ')
text(0.5,0.1,sprintf(' M = %0.2f',X(1,nSubjs+1)),'fontsize',20,'horizontalalignment','center')

%
% False Alarm Rate
axes('position',[0.23 0.1 0.15 0.85]); hold on;
xlim([0 1]); ylim([-0.1 1.1])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,X(2,ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD)
end
plot([0.2 0.8],[1 1]*X(2,nSubjs+1),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 0.5 1],'linewidth',2,'yticklabel',{0,'',1}); 
xlabel(' False Alarm Rate ')
text(0.5,0.4,sprintf(' M = %0.2f',X(2,nSubjs+1)),'fontsize',20,'horizontalalignment','center')

% d-Prime
axes('position',[0.41 0.1 0.15 0.85]); hold on;
xlim([0 1]); ylim([-0.4 4.4])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,X(3,ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD)
end
plot([0.2 0.8],[1 1]*X(3,nSubjs+1),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 2 4],'linewidth',2); 
xlabel(' d-prime ')
text(0.5,0.4,sprintf(' M = %0.2f',X(3,nSubjs+1)),'fontsize',20,'horizontalalignment','center')

% RT hits
axes('position',[0.61 0.1 0.15 0.85]); hold on;
xlim([0 1]); ylim([-0.2 2.2])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,X(5,ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD)
end
plot([0.2 0.8],[1 1]*X(5,nSubjs+1),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 1 2],'linewidth',2); 
xlabel(' RTs-Hits(s) ')
text(0.5,0.2,sprintf(' M = %0.2fs',X(5,nSubjs+1)),'fontsize',20,'horizontalalignment','center')

%
% RT CRs
axes('position',[0.79 0.1 0.15 0.85]); hold on;
xlim([0 1]); ylim([-0.2 2.2])

for ss= 1:nSubjs
    s=scatter(0.5+randn*0.1,X(6,ss),shapes(ss));
    set(s,'markerfaceColor','k','markeredgeColor','k','sizeData',sD)
end
plot([0.2 0.8],[1 1]*X(6,nSubjs+1),'linewidth',3,'color',0.4*ones(1,3))
set(gca,'fontsize',20,'xtick',[],'ytick',[0 1 2],'linewidth',2); 
xlabel(' RTs-CRs(s) ')
text(0.5,0.2,sprintf(' M = %0.2fs',X(6,nSubjs+1)),'fontsize',20,'horizontalalignment','center')

fN = 'retPerfScatter';
fP = '~/Google Drive/Research/ECoG_SS2e/';
cPath = pwd;
cd(fP)
addpath(cPath)
addpath([cPath '/Plotting/'])

print(gcf,'-dsvg',fN)
eval(['!' inkscapePath ' -z ' fN '.svg' ' --export-pdf=' fN '.pdf'])
eval(['!rm ' fN '.svg'])

cd(cPath)