

function parametricTimeCoursePlot(X,t,info)

inkscapePath='/Applications/Inkscape.app/Contents/Resources/bin/inkscape';


% defaults
if ~isfield(info,'yCenter')
    info.yCenter        = 0;
end
if ~isfield(info, 'smoother')
    info.smoother		= 'loess';
    info.smootherSpan 	= 0.15;
end

if ~isfield(info,'savePath')
    info.savePath       = pwd;
elseif ~exist(info.savePath,'dir')
    mkdir(info.savePath);
end

baseCols 	= [ 0.5 	0.0 	0.0;
				0.0 	0.0 	0.5;
				0.0 	0.5 	0.0;
				0.5 	0.5 	0.0;
				0.5 	0.0 	0.5]; 

endCols 	= [ 1.0 	0.8 	0.8;
				0.8 	0.8 	1.0;
				0.8 	1.0 	0.8;
				1.0 	1.0 	0.8;
				1.0 	0.8 	1.0];

nTypes 	= numel(X);
nLines 	= size(X{1},1);
cols 	= cell(nTypes,1);

if isfield(info,'colID')
    colIter = info.colID;
else
    colIter = 1:nTypes;
end

cnt = 1;
for jj = colIter
	cols{cnt} = [linspace(baseCols(jj,1),endCols(jj,1),nLines)' ...
				linspace(baseCols(jj,2),endCols(jj,2),nLines)' ...
				linspace(baseCols(jj,3),endCols(jj,3),nLines)'];
    cnt = cnt + 1;
end

if isfield(info,'yLimits')
    minY = info.yLimits(1);
    maxY = info.yLimits(2);
else
    maxY 	= max(cellfun(@max,X)); 
	maxY 	= maxY + 0.25*maxY;
	minY 	= min(cellfun(@min,X)); 
	minY 	= minY + 0.25*minY;
end

if ~isfield(info,'yRefLimits')
    info.yRefLimits = [minY*0.3 maxY*0.3];
end

figure(1); clf; hold on;
set(gcf,'position',[100 200 750 500],'paperpositionmode','auto')

ylim([minY maxY])
xlim([t(1)-0.01 t(end)+0.01])

if isfield(info,'toi') && isfield(info,'shade_toi')
    if info.shade_toi==1
        pp=patch([info.toi(1) info.toi(1) info.toi(2) info.toi(2)],...
            [minY maxY maxY minY],[0.92 0.92 0.92]);
        set(pp,'linestyle','none')
    end
end

plot(xlim,[info.yCenter info.yCenter],'--k','linewidth',2,'color',0.3*ones(3,1))
plot([0 0],info.yRefLimits,'--k','linewidth',2,'color',0.3*ones(3,1))

for jj = 1:nTypes
	for kk = 1:nLines
		if strcmp(info.smoother,'none')	
			x = X{jj}(kk,:);
		else
			x = smooth(X{jj}(kk,:),info.smootherSpan,info.smoother);		
		end

		plot(t,x,'color',cols{jj}(kk,:),'linewidth',2);
	end
end

set(gca,'linewidth',2)
set(gca,'fontweight','bold')
set(gca,'fontsize',15);
minY = minY + 0.05*abs(minY);
ylim([minY maxY])

x=get(gca,'YTick');
if (max(x)==-min(x))
    set(gca,'YTick',[min(x) 0  max(x)])
elseif (max(x) > -min(x)) && ( min(x) < 0)
    set(gca,'YTick',[min(x) 0 max(x)/2 max(x)])
elseif (max(x) > -min(x))
    set(gca,'YTick',[0 max(x)/2 max(x)])
elseif (min(x) < -max(x)) && ( max(x) > 0)
    set(gca,'YTick',[min(x) min(x)/2 0 max(x)])
elseif (min(x) < -max(x))
    set(gca,'YTick',[min(x) min(x)/2 0])
end
set(gca,'yTickLabel',get(gca,'yTick'));
hold off;

if isfield(info,'yAxisRightLoc')
	set(gca,'YAxisLocation','right')
end

if isfield(info,'xtick') && isfield(info,'xticklabel')
	set(gca,'xtick',info.xtick)
	set(gca,'xticklabel',info.xticklabel)
end

if isfield(info,'legend') & nTypes==2
    cols2=[cols{1};flipud(cols{2})];
    colormap(cols2)
    hh=colorbar;    
    set(hh,'box','off','ytick',[],'YColor',[1 1 1],'XColor',[1 1 1]);
end
set(gca,'fontSize',24)

cPath = pwd;
cd(info.savePath)
addpath(cPath)
addpath([cPath '/Plotting/'])

try
    fN = strcat(info.fileName,'.svg');
    plot2svg(fN,gcf)
    
    eval(['!' inkscapePath ' -z ' fN ' --export-pdf=' info.fileName '.pdf'])
    eval(['!rm ' fN])
catch
end
cd(cPath)