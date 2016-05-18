function f=renderChanWeights(elecLocs,Weights,CM)

% supplemental figure
% renders each brain in native and mni space
f = figure(); clf;
figW = 1000;
figH = 600;
set(gcf,'position',[100 200,figW,figH],'PaperPositionMode','auto','color','w')

resolution = 600;
%filename = 'MNIsChanRenderings';

%inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';
%load('~/Google Drive/Research/ECoG_SS2e/data_results/Renderings/electrodeLocs.mat')

view{1}      = [310,30];
view{2}      = [50,30];

lMNIcortex   =  elecLocs.lMNIcortex;
rMNIcortex   =  elecLocs.rMNIcortex;

%% render
clf;
ha = tight_subplot(1,2);

% left
axes(ha(1));
ctmr_gauss_plot(ha(1),lMNIcortex,[0 0 0],0,'l')
loc_view(view{1}(1),view{1}(2))
chans = find(elecLocs.hemChans==1);
el_add(elecLocs.MNILocs(chans,:),'k',40) 
for cc = 1:numel(chans)
    if ~isnan(Weights(chans(cc)))
        el_add(elecLocs.MNILocs(chans(cc),:),CM(Weights(chans(cc)),:),35);
    end
end

% right
axes(ha(2));
ctmr_gauss_plot(ha(2),rMNIcortex,[0 0 0],0,'r')
loc_view(view{2}(1),view{2}(2))
chans = find(elecLocs.hemChans==2);
el_add(elecLocs.MNILocs(chans,:),'k',40) 
for cc = 1:numel(chans)
    if ~isnan(Weights(chans(cc)))
        el_add(elecLocs.MNILocs(chans(cc),:),CM(Weights(chans(cc)),:),35);
    end
end
%print(gcf,'-dtiff',['-r' num2str(resolution)],[savePath filename])
