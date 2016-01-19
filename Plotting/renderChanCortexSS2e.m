function renderChanCortexSS2e()

% supplemental figure
% renders each brain in native and mni space
f = figure(1); clf;
figW = 1000;
figH = 600;
set(gcf,'position',[100 200,figW,figH],'PaperPositionMode','auto','color','w')

resolution = 600;
plotPath = '/Users/alexandergonzalez/Google Drive/Research/ECoG_SS2e/Plots/';
filename = 'MNIsChanRenderings';

%inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';
load('~/Google Drive/Research/ECoG_SS2e/data_results/electrodeLocs.mat')

view{1}      = [310,30];
view{2}      = [50,30];

%hemId       = {'l','r'};

lMNIcortex   =  elecLocs.lMNIcortex;
rMNIcortex   =  elecLocs.rMNIcortex;

% ROI channel colors
colors      = [];
colors{1}  = [0.9 0.2 0.2];
colors{2}  = [0.1 0.5 0.8];
colors{3}   = [0.2 0.6 0.3];

%% render
clf;
ha = tight_subplot(1,2);

axes(ha(1));
ctmr_gauss_plot(ha(1),lMNIcortex,[0 0 0],0,'l')
loc_view(view{1}(1),view{1}(2))
el_add(elecLocs.MNILocs(elecLocs.hemChans==1,:),'k',40) 
for rr = 1:3
    el_add(elecLocs.MNILocs(elecLocs.ROIid==rr&elecLocs.hemChans==1,:),colors{rr},35);
end

axes(ha(2));
ctmr_gauss_plot(ha(2),rMNIcortex,[0 0 0],0,'r')
loc_view(view{2}(1),view{2}(2))
el_add(elecLocs.MNILocs(elecLocs.hemChans==2,:),'k',40) 
for rr = 1:3
    el_add(elecLocs.MNILocs(elecLocs.ROIid==rr&elecLocs.hemChans==2,:),colors{rr},35);
end

 print(gcf,'-dtiff',['-r' num2str(resolution)],[plotPath filename])
