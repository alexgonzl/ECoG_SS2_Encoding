
col = 'b'

for ss = 1:5
seed = ss;
rng(seed);

nPatterns = [6 6];

switch col 
    case 'g'
        mainCol     = [0.4 0.7 0.4];
    case 'r'
        mainCol     = [0.8 0.4 0.4];
    case 'o'
        mainCol     = [150 120 80]/255;
    case 'c'
        mainCol     = [100 180 180]/255;
    case 'gray'
        mainCol     = [200 200 200]/255;
    case 'b'
        mainCol     = [100 150 200]/255;
end

figure(1); clf;
set(gcf,'paperpositionmode','auto','color','none')
set(gcf,'paperUnits','points','papersize',[400 400],'paperposition',[0 0 400 400])
set(gcf,'position',[100,100,400,400]); hold on;


for ii = 1:nPatterns(1)
    for jj = 1:nPatterns(2)
        s = scatter(ii,jj,'o');
        s.SizeData = 3e3;
        s.MarkerFaceAlpha=1;
        s.MarkerFaceColor = min(mainCol + randn(1,3)*0.1,1);
        s.MarkerEdgeColor='none';
    end
end
axis off
axis square
set(gca,'color','none')

fPath = '/Users/alexandergonzalez/Google Drive/Research/Thesis/Defense/figures/';
print(gcf,'-dpdf', [fPath 'Patterns_' col num2str(seed)])
end