
X = permute(data.BinTStatAvg,[2 1 3]);
Y = X(:,:);

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(Y,'econ');
%%
MSE = size(COEFF,2);
for jj = 1:size(COEFF,2)
    MSE(jj)=sum(sum(abs(SCORE(:,1:jj)*COEFF(:,1:jj)'-Y).^2))./numel(Y);
end


%%
GMModels = cell(10,1); % Preallocation
options = statset('MaxIter',1000);
rng(1); % For reproducibility

for j = 1:10
    GMModels{j} = fitgmdist(SCORE,j,'Options',options,'Replicates',200,'RegularizationValue',0.4);
    fprintf('\n GM Mean for %i Component(s)\n',j)
    Mu = GMModels{j}.mu
    BIC(j)= GMModels{j}.BIC;
    AIC(j)= GMModels{j}.AIC;
end
[minBIC,numComponentsB] = min(BIC);
[minAIC,numComponentsA] = min(AIC);
numComponentsA
numComponentsB
%%x
figure
for j = 1:3
    subplot(2,2,j)
    gscatter(SCORE(:,1),SCORE(:,2),data.ROIid)
    h = gca;
    hold on
    ezcontour(@(x1,x2)pdf(GMModels{j},[x1 x2]),...
        [h.XLim h.YLim],100)
    title(sprintf('GM Model - %i Component(s)',j));
    xlabel('1st principal component');
    ylabel('2nd principal component');
    if(j ~= 3)
        legend off;
    end
    hold off
end
g = legend;
g.Position = [0.7 0.25 0.1 0.1];