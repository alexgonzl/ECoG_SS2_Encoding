function ha = subSpecPlot(dataStruct)

%% set colormap
colMapFun = @(bot,top,N)([ linspace(top(1),bot(1),N)', ...
    linspace(top(2),bot(2),N)', linspace(top(3),bot(3),N)']);

col1 = [0.2081    0.1663    0.5292];%dataStruct.col1;%[0.8 0.5 0.2];
col2 = [1 1 1];
col3 = [0.9763    0.9831    0.0538];%dataStruct.col2;%[0.2 0.4 0.8];
nLevels = 50;

colMap = [colMapFun(col2,col3,nLevels);
                col2;
                colMapFun(col1,col2,nLevels)];

%% set axes location
tcW = 0.8;
tcH = 0.8;

ha(1)=axes('position',[0.1  0.1 tcW tcH]);

%% parameters
freqs 	= dataStruct.freqs; 
nFreqs 	= numel(freqs);

time 	= dataStruct.time;
X 		= dataStruct.data;
Zlimits = dataStruct.Zlims;
timeTicks = {-0.2:0.2:1};

X(~(abs(X)>2))=nan;
% first spectrogram
axes(ha(1))
imagesc(time,1:nFreqs,X,Zlimits);
hold on;
plot([0 0],ylim,'--k','linewidth',2)

set(gca,'ytick',1:nFreqs)
set(gca,'yticklabel',freqs(get(gca,'ytick')))
set(gca,'fontsize',14,'fontWeight','normal')
set(gca,'xtick',timeTicks{1},'xticklabel',[])
set(gca,'lineWidth',0.5)
axis xy
colormap(gca,colMap)

if isfield(dataStruct,'ylabel')
	ylabel(dataStruct.ylabel)
end

colormap(gca,colMap)