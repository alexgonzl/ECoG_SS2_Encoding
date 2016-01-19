function renderChansSS2e(chansIdx,AnalysisType)

% supplemental figure
% renders each brain in native and mni space
f = figure(1); clf;
figW = 1000;
figH = 600;
set(gcf,'position',[100 200,figW,figH],'PaperPositionMode','auto','color','w')

resolution = 600;
plotPath = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Plots/';
filename = strcat('MNIsChanRenderings',AnalysisType );

%inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';
load('~/Google Drive/Research/ECoG_SS2e/data_results/electrodeLocs.mat')

view{1}      = [310,30];
view{2}      = [50,30];

%hemId       = {'l','r'};

lMNIcortex   =  elecLocs.lMNIcortex;
rMNIcortex   =  elecLocs.rMNIcortex;

% ROI channel colors
colors      = [];
colors{1}  = [0.5 0.6 0.7];
colors{2}  = [0.7 0.6 0.5];
colors{3}  = [0.6 0.7 0.6];

%% render
clf;
ha = tight_subplot(1,2);

axes(ha(1));
ctmr_gauss_plot(ha(1),lMNIcortex,[0 0 0],0,'l')
loc_view(view{1}(1),view{1}(2))
el_add(elecLocs.MNILocs(elecLocs.hemChans==1,:),'k',40) 
for rr = 1:3
    el_add(elecLocs.MNILocs(chansIdx==rr&elecLocs.hemChans==1,:),colors{rr},35);
end

axes(ha(2));
ctmr_gauss_plot(ha(2),rMNIcortex,[0 0 0],0,'r')
loc_view(view{2}(1),view{2}(2))
el_add(elecLocs.MNILocs(elecLocs.hemChans==2,:),'k',40) 
for rr = 1:3
    el_add(elecLocs.MNILocs(chansIdx==rr&elecLocs.hemChans==2,:),colors{rr},35);
end
% set legend axes
axes('position',[0.94, 0.45, 0.03 0.1]); hold on;
plot([0.2],[3],'o','color',colors{1},'markersize',15,'markerfacecolor',colors{1})
plot([0.2],[2],'o','color',colors{2},'markersize',15,'markerfacecolor',colors{2})
plot([0.2],[1],'o','color',colors{3},'markersize',15,'markerfacecolor',colors{3})

xlim([0 0.8]); ylim([0.5 3.5])
text(0.5,3,'CL1','fontsize',16)
text(0.5,2,'CL2','fontsize',16)
text(0.5 ,1,'CL3','fontsize',16)
set(gca,'visible','off')

 print(gcf,'-dtiff',['-r' num2str(resolution)],[plotPath filename])
