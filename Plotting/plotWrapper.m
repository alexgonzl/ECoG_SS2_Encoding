function plotWrapper(X,t,info)

inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';


% defaults
if ~isfield(info,'centerType')
	info.centerType 	= 'mean';
end
if ~isfield(info,'yLimits')
	info.yLimits        = [];
end
if ~isfield(info,'yCenter')
	info.yCenter        = 0;
end
if ~isfield(info, 'smoother')
	info.smoother		= 'loess';
	info.smootherSpan 	= 0.10;
end	
if ~isfield(info,'errorBars')
	info.errorBars      = 'se';
end
if ~isfield(info,'cols')
	info.cols           = 'rbgk';
end
if ~isfield(info,'savePath')
	info.savePath       = pwd;
elseif ~exist(info.savePath,'dir')
	mkdir(info.savePath);
end

figure(1); clf;
set(gcf,'position',[100 200 750 500],'paperpositionmode','auto')
h = plotNTraces(X,t,'cols',info.cols,...
		'smoother',info.smoother,'smootherSpan',info.smootherSpan,...
		'centerType',info.centerType,'errorBars',info.errorBars,...
		'yLimits',info.yLimits);

 if isfield(info,'legend')
     l = legend([h.h1.mainLine h.h2.mainLine], info.legend,'location','best');
     set(l,'box','off')
 end
 
cPath = pwd;
cd(info.savePath)
addpath(cPath)
addpath([cPath '/Plotting/'])

try
fN = strcat(info.fileName,'.svg');
plot2svg(fN,h.f)

eval(['!' inkscapePath ' -z ' fN ' --export-pdf=' info.fileName '.pdf'])
eval(['!rm ' fN])
catch
end
cd(cPath)