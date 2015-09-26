
function h = plotNTraces(X,t,varargin)
% h = plotNTraces(X,t,colors,smoother,smootherSpan)
% function that plots traces with mean and cells
% takes a cell X: number of cells indicate the number of traces; each
% cell should be a matrix with:
% 1st dim = observations
% 2nd dim = time
% t = time units
% colors for each trace
% smoother: {'loess','lowess','moving'}
% smootherSpan: 0-1 for loess, nsamps for moving
% centerType indicates if it should use the mean or the median


% defaults
cols          = 'ocfbrkly';
centerType      = 'mean';
errorBars       = 'se';

% optional inputs
nOptInputs  = numel(varargin);
if nOptInputs >0
    assert(mod(nOptInputs,2)==0,'unpaired input')
    assert(all(ismember(varargin(1:2:end)',...
        {'cols','smoother','smootherSpan','centerType',...
        'yLimits','errorBars','yCenter','yRefLimits'})),...
         'invalid input')
    
    for oo=1:2:nOptInputs
        if ~isempty(varargin{oo+1})
            eval([varargin{oo} '= varargin{oo+1};'])
        end
    end
end

N = numel(X);
nSamps = size(X{1},2);

% colors 
o = [0.9 0.5 0.1]; % orange
c = [0.1 0.8 0.9]; % cyan
g = [0.2 0.6 0.3]; % green
b = [0.1 0.5 0.8]; % blue
r = [0.9 0.2 0.2]; % red
k = [0.05 0.05 0.05]; % black
l = [0.5 0.2 0.5]; % magenta
y = [0.9 0.9 0.3]; % yellow

colors = zeros(N,3);

n   = zeros(N,1);
m   = zeros(N,nSamps);
eb  = zeros(N,nSamps);

for i = 1:N
    colors(i,:) = eval(cols(i));
    x = X{i};    
    m(i,:)  = eval([centerType '(x)']);
    n(i)    = size(x,1);
    eb(i,:) = eval([errorBars '(x)']);
    if exist('smoother','var') 
        if ~strcmp(smoother,'none')
            m(i,:)  = smooth(m(i,:),smootherSpan,smoother);
            eb(i,:) = smooth(eb(i,:),smootherSpan,smoother);
        end
    end
end

h.f=figure(gcf);cla;hold on;

if exist('yLimits','var')
    minY = yLimits(1);
    maxY = yLimits(2);
else
    minY = min(min(m-eb));
    minY = minY - 0.25*abs(minY);
    maxY = max(max(m+eb));
    maxY = maxY + 0.25*abs(maxY);    
end

if exist('yRefLimits','var')
    if isempty(yRefLimits)
        yRefLimits = [minY*0.3 maxY*0.3];
    end
else
    yRefLimits = [minY*0.3 maxY*0.3];
end

if ~exist('yCenter','var')
    yCenter = 0;
end
    
ylim([minY maxY])

for c = 1:N
    plot(t,m(c,:)+eb(c,:),'color',colors(c,:));
    plot(t,m(c,:)-eb(c,:),'color',colors(c,:));
end

xlim([t(1)-0.01 t(end)+0.01])
plot(xlim,[yCenter yCenter],'--k','linewidth',2,'color',0.3*ones(3,1))
plot([0 0],yRefLimits,'--k','linewidth',2,'color',0.3*ones(3,1))

for c = 1:N
    h.(['h' num2str(c)]) = shadedErrorBar(t,m(c,:),eb(c,:),{'color', colors(c,:),'linewidth',3},1);
    set(h.(['h' num2str(c)]).mainLine,'marker','.','MarkerSize',11)
end
set(gca,'linewidth',2)
set(gca,'fontweight','bold')
set(gca,'fontsize',16);
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
