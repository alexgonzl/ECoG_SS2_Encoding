
Ks = 20;
AIC = zeros(Ks,1);
BIC = zeros(Ks,1);
%X=permute(out.Projections,[2 1 3]); X=X(:,:)';
X = Y;
N=size(X,1);
P=size(X,2);
for ii = 1:Ks
    m = fitgmdist(X,ii,'RegularizationValue',1/P,'replicates',10);
    AIC(ii) = m.AIC;
    BIC(ii) = m.BIC;
    fprintf('GMM Fit completed for %i  Nodes \n',ii)
end
%% ii=10 is the best solution

m = fitgmdist(X,3,'RegularizationValue',1/P,'replicates',10);