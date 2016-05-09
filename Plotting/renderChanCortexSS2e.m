function renderChanCortexSS2e(savePath)

% supplemental figure
% renders each brain in native and mni space
f = figure(1); clf;
figW = 1000;
figH = 600;
set(gcf,'position',[100 200,figW,figH],'PaperPositionMode','auto','color','w')

resolution = 600;
filename = 'MNIsChanRenderings';

%inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';
load('~/Google Drive/Research/ECoG_SS2e/data_results/Renderings/electrodeLocs.mat')

view{1}      = [310,30];
view{2}      = [50,30];

%hemId       = {'l','r'};

lMNIcortex   =  elecLocs.lMNIcortex;
rMNIcortex   =  elecLocs.rMNIcortex;

% ROI channel colors
nROIs       = 5;
colors      = [];
colors{1}  = [240 35 17]/255;
colors{2}  = [2 93 140]/255;
colors{3}  = [122 179 23]/255;
colors{4}  = [0.7 0.4 0.8];
colors{5}  = [0.7 0.7 0.4];

%% render
clf;
ha = tight_subplot(1,2);

axes(ha(1));
ctmr_gauss_plot(ha(1),lMNIcortex,[0 0 0],0,'l')
loc_view(view{1}(1),view{1}(2))
el_add(elecLocs.MNILocs(elecLocs.hemChans==1,:),'k',40) 
for rr = 1:nROIs
    el_add(elecLocs.MNILocs(elecLocs.ROIid==rr&elecLocs.hemChans==1,:),colors{rr},35);
end

axes(ha(2));
ctmr_gauss_plot(ha(2),rMNIcortex,[0 0 0],0,'r')
loc_view(view{2}(1),view{2}(2))
el_add(elecLocs.MNILocs(elecLocs.hemChans==2,:),'k',40) 
for rr = 1:nROIs
    el_add(elecLocs.MNILocs(elecLocs.ROIid==rr&elecLocs.hemChans==2,:),colors{rr},35);
end

 print(gcf,'-dtiff',['-r' num2str(resolution)],[savePath filename])
