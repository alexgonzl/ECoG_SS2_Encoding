function out = getGroupStats(X,Y)
% get ttest statistiscts for X, for Y, and X against Y
% rows -> observations 
% columns -> variables
% X and Y must have the same number of columns.


P = size(X,2);
out = zeros(3,P,2);

[~,p,~,t]=ttest(X);
out(1,:,1) = p;
out(1,:,2) = t.tstat;

[~,p,~,t]=ttest(Y);
out(2,:,1) = p;
out(2,:,2) = t.tstat;

[~,p,~,t]=ttest2(X,Y);
out(3,:,1) = p;
out(3,:,2) = t.tstat;