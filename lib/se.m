function out = se(x)

n = size(x,1);
x(isnan(x))=[];
out = std(x)/sqrt(n-1);